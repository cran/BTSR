subroutine linkR(link, par, n, ind, x, gx, dlink)
  !------------------------------------------------------------------------
  !  Interface for R
  !  Evaluates:
  !     - g(y) = the link function at y, if id(1) = 1
  !     - g^{-1}(gy) = the inverse of the link function at gy, if id(2) = 1
  !     - dg(y)/dy = the derivative of the link function at y, if id(3) = 1
  !
  !  Input:
  !   link : an integer indicating the type of link
  !             0 = polynomial --> g(y) = a * y**b
  !             1 = logit      --> g(y) = log((y - ymin) / (ymax - y))
  !             2 = logarithm  --> g(y) = log(y)
  !             3 = loglog     --> g(y) = log(-log((y - ymin) / (ymax - ymin)))
  !             4 = c-loglog   --> g(y) = log(-log(1.d0 - (y - ymin) / (ymax - ymin)))
  !             5 = SIP        --> g(y) = 1 / (a + y)**b, a = 0 or 1
  !   par  : vector of size 4 containing
  !           - the lower and upper limits for y
  !           - link parameters a and b.
  !   n    : sample size
  !   ind  : the indicator of which quantity to calculate
  !   x    : the values of the variable y
  !   gx   : the values of g(y)
  !   dlink: the values of dg(y)/dy
  !------------------------------------------------------------------------
  ! Last revision: February, 2025
  !   - replaced the size 2 vector 'ylim' by the size 4 vector 'par'
  use base
  implicit none
  integer :: n, link, ind(3)
  real(dp) :: x(n), par(4)
  real(dp) :: gx(max(1, n * (1 - ind(3))))
  real(dp) :: dlink(max(1, n * ind(3)))
  type(argslink) :: args
  integer :: i
  args%link = link
  args%lower = par(1)
  args%upper = par(2)
  args%par = par(3:4)

  if(ind(1) == 1) then
     do i = 1, n
        gx(i) = linkfun(x(i), args)
     end do
     return
  else if (ind(2) == 1) then
     do i = 1, n
        x(i) = linkinv(gx(i), args)
     end do
     return
  else
     do i = 1,n
        dlink(i) = diflink(x(i), args)
     end do
     return
  end if
end subroutine linkR

subroutine simbtsr(code, length, order, ts, xreg1, xreg2, tstart, xstart, link, lconfig,  &
     np, pdist, xregar, alpha, beta, phi, theta, d, rev)
  use base
  !---------------------------------------------------------------------
  !  Interface for R
  !  Simulates a BTSR model
  !
  ! Input:
  !  code    : the model's code
  !            (see the list of current models in the main file)
  !  length  : a vector of size 2 with burn in and final sample size
  !  order   : a 2 by 4 matrix with (nreg, p, q, inf)
  !  ts      : the simulated time series yt, mut, eta1t, e1t, nut, gnut, eta2t
  !  xreg1   : matrix of regressors corresponding to part 1
  !  xreg2   : matrix of regressors corresponding to part 2
  !  tstart  : starting value for y, g(nu) and g23(e1t).
  !  xstart  : matrix of maxval(nreg) rows by 2 rows with starting values for the regressors
  !            (to be used when t < 1)
  !  link    : a vector of size 8, indicating the links to be used in the model
  !  lconfig : a matrix of size 8 by 4 with
  !             - the lower and upper bounds related to the link functions
  !             - the constant and power arguments for polynomial/SIP links
  !  np      : number of parameters in pdist
  !  pdist   : nuisance parameters to the conditional distribution
  !  xregar  : indicates if xreg must be included in the AR recursion
  !  alpha   : a vector of size 2 with the intercept for each part
  !  beta    : a matrix with 2 rows (regressors' coefficients)
  !  phi     : a matrix with 2 rows (AR coefficients)
  !  theta   : a matrix with 2 rows (MA coefficients)
  !  d       : a vector of size 2 with the long memory parameter
  !  rev     : flag to return error code
  !---------------------------------------------------------------------
  ! Implemented by Taiane S. Prass - PPGEst / UFRGS
  ! February, 2024 - Porto Alegre
  !
  ! Last update: February, 2025
  !  - simplified the imput using matrices for order and the time series
  !  - added the rng type to eliminate interfaces
  !---------------------------------------------------------------------
  implicit none
  integer  :: code
  integer, target  :: length(2), np
  integer, pointer :: n, burn
  integer  :: order(2,4)
  integer  :: link(8), xregar(2)
  integer  :: rev
  real(dp) :: ts(sum(length), 8)
  real(dp) :: pdist(np), alpha(2), d(2)
  real(dp) :: lconfig(8, 4), tstart(3)
  real(dp) :: xstart(2, maxval([1, order(:,1)]))
  real(dp) :: beta(2, maxval([1, order(:,1)]))
  real(dp) :: phi(2, maxval([1, order(:,2)]))
  real(dp) :: theta(2, maxval([1, order(:,3)]))
  real(dp) :: xreg1(sum(length), max(1, order(1,1)))
  real(dp) :: xreg2(sum(length), max(1, order(2,1)))
  type(argsDist) :: rdist

  burn => length(1)
  n    => length(2)

  ! Initialize the argsModel object
  call rdist%init(current_models(code))  ! Set the model

  call sim_model(rdist, n, burn, np, pdist, onlymu(code), order,  xregar, alpha, beta, &
       phi, theta, d, link, lconfig, xreg1, xreg2, xstart, tstart, ts, rev)
  return
end subroutine simbtsr

subroutine extractbtsr(code, length, order, ts, xreg1, xreg2, tstart, xstart, xnew1, xnew2,  &
     forecast, link, lconfig, npar, par, xregar, nfix, alpha, flagsb, beta, flagsar, &
     phi, flagsma, theta, d, np, pdist, extras, sll, U, K, nd, Dg, T, E, h, ierr)
  use base
  !-----------------------------------------------------------------------------
  !
  !  Extracting components from a BTSR model.
  !  Wrapper to call the subroutines from R.
  !
  ! Input:
  !  code    : the model's code
  !            (see the list of current models in the main file)
  !  length  : size 2 vector with n = sample size and nnew =  number of forecasts
  !  order   : a 2 by 4 matrix with (nreg, p, q, inf)
  !  ts      : the observed time series yt and the extracted time series
  !            mut, eta1t, e1t, nut, gnut, eta2t, e2t
  !  xreg1   : matrix of regressors corresponding to part 1
  !  xreg2   : matrix of regressors corresponding to part 2
  !  tstart  : starting value for y, g2(nu) and g23(e1). To be used when t < 1
  !  xstart  : matrix of maxval(nreg) rows by 2 rows with starting values for the regressors
  !            (to be used when t < 1)
  !  xnew1   : matrix with new values of regressors corresponding to part 1
  !  xnew2   : matrix with new values of regressors corresponding to part 2
  !  forecast: matrix with the forecast for mut, nut and varthetat
  !  link    : a vector of size 8, indicating the links to be used in the model
  !  lconfig : a matrix of size 8 by 4 with
  !             - the lower and upper bounds related to the link functions
  !             - the constant and power arguments for polynomial/SIP links
  !  npar   : number of non-fixed parameters for each part
  !  coefs  : non-fixed parameters
  !  xregar : indicates if xreg must be included in the AR recursion
  !  nfix   : a 2 by 5 matrix with the number of fixed parameters
  !  alpha  : a vector of size 2 with the intercept for each part
  !  flagsb : the lags of beta that must be fixed
  !  beta   : a matrix with 2 rows (regressors' coefficients)
  !  flagsar: the lags of phi that must be fixed
  !  phi    : a matrix with 2 rows (AR coefficients)
  !  flagsma: the lags of theta that must be fixed
  !  theta  : a matrix with 2 rows (MA coefficients)
  !  d      : a vector of size 2 with the long memory parameter
  !  np     : number of parameters in pdist
  !  pdist  : nuisance parameters to the conditional distribution
  !  extras  : a vector containing
  !           - m: an integer indicating where the sum of the log-likelihood starts
  !           - llk: indicates if the sum of the log-likelihood must be returned
  !           - sco: indicates if the score vector must be returned
  !           - info: indicates if the information matrix must be returned
  !           - extra: indicates if the extra matrices and vectors must be returned
  !  sll    : the sum of the log-likelihood
  !  U      : the score vector
  !  K      : the information matrix
  !  nd     : number of columns in D.
  !  D      : matrix with derivatives deta1/drho, deta2/drho and deta2/dlambda
  !  T      : 2 by b matrix with derivatives dmu/deta1 and dnu/deta2
  !  E      : 3 by n matrix with the expected values E(dl2/dmu2), E(dl2/dmUdnu) and E(dl2/dnu2)
  !  h      : 2 by n matrix with derivatives dl/dmu, dl/dnu
  !  ierr   : error code
  !---------------------------------------------------------------------
  ! Implemented by Taiane S. Prass - PPGEst / UFRGS
  ! April, 2024 - Porto Alegre
  !
  ! Last update: March, 2025
  !  - simplified the imput using matrices for order and the time series
  !  - added the rng type to eliminate interfaces
  !  - added ierr
  !---------------------------------------------------------------------
  implicit none
  ! time series
  integer  :: code, nd, ierr
  integer  :: length(2), link(8), order(2, 4)
  real(dp) :: lconfig(8, 4)
  real(dp) :: ts(length(1), 10)
  real(dp) :: xreg1(length(1), max(1, order(1, 1)))
  real(dp) :: xreg2(length(1), max(1, order(2, 1)))
  real(dp) :: tstart(3), xstart(2, maxval([1, order(:, 1)]))
  real(dp) :: xnew1(max(1, length(2)), max(1, order(1, 1)))
  real(dp) :: xnew2(max(1, length(2)), max(1, order(2, 1)))
  real(dp) :: forecast(max(1, length(2)), 5)
  ! parameters
  integer  :: np, npar(2), xregar(2), nfix(2, 5)
  integer  :: flagsb(2, maxval([1, nfix(:, 2)]))
  integer  :: flagsar(2, maxval([1, nfix(:, 3)]))
  integer  :: flagsma(2, maxval([1, nfix(:, 4)]))
  real(dp) :: par(sum(npar)), alpha(2), d(2), pdist(np)
  real(dp) :: beta(2, maxval([1, nfix(:, 2)]))
  real(dp) :: phi(2, maxval([1, nfix(:, 3)]))
  real(dp) :: theta(2, maxval([1, nfix(:, 4)]))
  ! likelihood
  integer  :: extras(5)
  real(dp) :: sll, U(max(1, sum(npar) * extras(3)))
  real(dp) :: K(max(1, sum(npar) * extras(4)), max(1, sum(npar) * extras(4)))
  real(dp) :: T(max(length(1) * extras(5), 1), max(2 * extras(5), 1))
  real(dp) :: h(max(length(1) * extras(5), 1), max(2 * extras(5), 1))
  real(dp) :: Dg(max(length(1) * extras(5), 1), max(nd * extras(5), 1))
  real(dp) :: E(max(length(1) * extras(5), 1), max(3 * extras(5), 1))
  type(argsModel) :: model

  ! Initialize the argsModel object
  call model%argsD%init(current_models(code))  ! Set the model

  model%argsD%lower = pdist(np-1) ! a
  model%argsD%upper = pdist(np)   ! b

  select case (trim(adjustl(current_models(code))))
  case("kuma")
     model%argsD%par = pdist(1) ! rho
  case("uw")
     model%argsD%par = pdist(1) ! rho
  case default
     model%argsD%par  = 0.d0
  end select

  ! allocating matrices and vectors and setting variables fixed values
  call get_model(model, length(1), order, ts(:,1), ts(:,2:3), xreg1, xreg2, tstart(1:2),  &
       xstart, link, lconfig, 0, npar, par, xregar, nfix, alpha, flagsb, beta, flagsar, phi, &
       flagsma, theta, d, extras(1:4), .true., ierr)

  ! Summary
  call final_model(model, sum(npar), par, length, ts(:,4:10), order(:,1), xnew1, xnew2, &
       forecast, order(:,4), sll, extras(3:5), U, K, nd, Dg, T, E, h)
  return
end subroutine extractbtsr

subroutine optimbtsr(method, code, length, order, ts, xreg1, xreg2, tstart, xstart, xnew1, xnew2, &
     forecast, link, lconfig, npar, par, nbd, bounds, xregar, nfix, alpha, flagsb, beta, flagsar, &
     phi, flagsma, theta, d, np, pdist, extras, sll, U, K, nd, Dg, T, E, h, cf1, nc2, cf2, neval, conv)
  use base
  implicit none
  !--------------------------------------------
  ! iprint < 0 = no print
  ! stopcr = stopping critereon  (1.d-4)
  !
  ! conv = 0 FOR SUCCESSFUL TERMINATION
  !      = 1 IF MAXIMUM NO. OF FUNCTION EVALUATIONS EXCEEDED
  !      = 2 IF NOP < 1 .or. STOPCR <= 0
  !---------------------------------------------
  ! time series
  integer  :: method, code, nd
  integer  :: length(2), link(8), order(2, 4)
  real(dp) :: lconfig(8, 4)
  real(dp) :: ts(length(1), 10)
  real(dp) :: xreg1(length(1), max(1, order(1, 1)))
  real(dp) :: xreg2(length(1), max(1, order(2, 1)))
  real(dp) :: tstart(3), xstart(2, maxval([1, order(:,1)]))
  real(dp) :: xnew1(max(1, length(2)), max(1, order(1, 1)))
  real(dp) :: xnew2(max(1, length(2)), max(1, order(2, 1)))
  real(dp) :: forecast(max(1, length(2)), 5)
  ! parameters
  integer  :: np, npar(2), nfix(2,5), xregar(2)
  integer  :: flagsb(2, maxval([1, nfix(:, 2)]))
  integer  :: flagsar(2, maxval([1, nfix(:, 3)]))
  integer  :: flagsma(2, maxval([1, nfix(:, 4)]))
  real(dp) :: par(sum(npar)), alpha(2), d(2), pdist(np)
  real(dp) :: beta(2, maxval([1, nfix(:, 2)]))
  real(dp) :: phi(2, maxval([1, nfix(:, 3)]))
  real(dp) :: theta(2, maxval([1, nfix(:, 4)]))
  ! optimization
  integer  :: extras(5), nbd(sum(npar))
  integer  :: cf1(2), neval, conv, nc2
  real(dp) :: bounds(sum(npar), 2), cf2(nc2)
  real(dp) :: sll, U(max(1, sum(npar) * extras(3)))
  real(dp) :: K(max(1, sum(npar) * extras(4)), max(1, sum(npar) * extras(4)))
  real(dp) :: T(max(length(1) * extras(5), 1), max(2 * extras(5), 1))
  real(dp) :: h(max(length(1) * extras(5), 1), max(2 * extras(5), 1))
  real(dp) :: Dg(max(length(1) * extras(5), 1), max(nd * extras(5), 1))
  real(dp) :: E(max(length(1) * extras(5), 1), max(3 * extras(5), 1))
  type(argsModel) :: model
  type(optimFunc) :: loglik
  type(optimizer) :: optim
  integer :: extr(4), cf(3)

  ! Initialize the argsModel object
  call model%argsD%init(current_models(code))  ! Set the model

  model%argsD%lower = pdist(np-1) ! a
  model%argsD%upper = pdist(np)   ! b

  select case (trim(adjustl(current_models(code))))
  case("kuma")
     model%argsD%par = pdist(1) ! rho
  case("uw")
     model%argsD%par = pdist(1) ! rho
  case default
     model%argsD%par  = 0.d0
  end select

  select case (method)
  case(0)
     ! L-BFGS-B
     extr = [extras(1), 1, 1, 0]  ! m, llk, sco, info
     cf = [cf1, 1] ! cf(3) is a multiplier for the size of the score vector
     loglik%loglik => loglik_dist
     optim%optim   => optim_lbfgsb
  case(1)
     ! Nelder-Mead
     extr = [extras(1), 1, 0, 0]  ! m, llk, sco, info
     cf = [cf1, 0] ! cf(3) is a multiplier for the size of the score vector
     loglik%functn => loglik_dist_nelder
     optim%optim   => optim_nelder
  end select

  ! allocating matrices and vectors and setting variables fixed values
  ! score vector related matrices will not be allocated at this point
  call get_model(model, length(1), order, ts(:,1), ts(:,2:3), xreg1, xreg2, tstart(1:2),  &
       xstart, link, lconfig, 0, npar, par, xregar, nfix, alpha, flagsb, beta, flagsar, phi, &
       flagsma, theta, d, extr, .true., conv)
  if(conv > 0) return

  ! Optimization subroutine
  call optim%optim(loglik, model, sum(npar), par, nbd, bounds, sll, U, cf, nc2, cf2, neval, conv)

  ! Summary of the final model
  if (extras(3) + extras(4) > 0) call allocate_SI(model, model%SI)
  call final_model(model, sum(npar), par, length, ts(:,4:10), order(:,1), xnew1, xnew2, &
       forecast, order(:,4), sll, extras(3:5), U, K, nd, Dg, T, E, h)
  return
end subroutine optimbtsr

subroutine predictbtsr(length, order, ts, xreg1, xreg2, xnew1, xnew2, forecast, &
     link, lconfig, npar, par,  xregar, nfix, alpha, flagsb, fvbeta, flagsar, fvar, &
     flagsma, fvma, d)
  !-----------------------------------------------------------------------------------------
  !
  !  Subroutine used for prediction.
  !  The values of y, g11y, g12y, eta, error,
  !
  !-----------------------------------------------------------------------------------------
  ! July, 2023
  !   - changed d == 0.0d0 to abs(d) < epsilon(1.d0)
  ! February, 2024
  !   - removed nu and fixnu
  !   - removed ylower and yupper and added lconfig
  ! Last revision: February, 2025
  !   - renamed gy to g12y
  use base
  implicit none
  integer  :: length(2), link(8), order(2,4), nfix(2,5)
  integer  :: npar(2), xregar(2)
  integer  :: flagsb(2, maxval([1, nfix(:, 2)]))
  integer  :: flagsar(2, maxval([1, nfix(:, 3)]))
  integer  :: flagsma(2, maxval([1, nfix(:, 4)]))
  real(dp) :: ts(length(1), 8)
  real(dp) :: xreg1(length(1), max(1, order(1, 1)))
  real(dp) :: xreg2(length(1), max(1, order(2, 1)))
  real(dp) :: xnew1(max(1, length(2)), max(1, order(1, 1)))
  real(dp) :: xnew2(max(1, length(2)), max(1, order(2, 1)))
  real(dp) :: forecast(max(1, length(2)), 5)
  real(dp) :: lconfig(8, 4)
  real(dp) :: par(sum(npar)), alpha(2), d(2)
  real(dp) :: fvbeta(2, maxval([1, nfix(:, 2)]))
  real(dp) :: fvar(2, maxval([1, nfix(:, 3)]))
  real(dp) :: fvma(2, maxval([1, nfix(:, 4)]))
  type(argsModel) :: model
  real(dp) :: xstart(2, maxval([1, order(:,1)]))
  integer  :: i, ierr

  !  If d = 0 uses inf = q
  do i = 1, 2
     !  If d = 0 uses inf = q
     if (abs(d(i)) < epsmch) then
        order(i, 4) = order(i, 3)
     else
        order(i, 4) = max(order(i, 4), order(i, 3))
     end if
  end do

  xstart = 0.d0     ! dummy

  ! allocating matrices and vectors and setting variables fixed values
  call get_model(model, length(1), order, ts(:,1), ts(:,2:3), xreg1, xreg2, [0.d0, 0.d0],  &
       xstart, link, lconfig, 0, npar, par, xregar, nfix, alpha, flagsb, fvbeta, flagsar, fvar, &
       flagsma, fvma, d, [0,0,0,0], .false., ierr)

  ! setting the values of eta1 and error (from 1 to n)
  model%cts(1)%et = ts(:,4)   ! e1t
  model%cts(1)%eta = 0.d0     ! g11(mut)
  model%cts(2)%estart = 0.d0  ! dummy
  model%cts(2)%gw =  ts(:,5)  ! gnut
  model%cts(2)%eta = ts(:,6)  ! g12(gnut)
  model%cts(2)%gi2 = ts(:,7)  ! g22(gnut)
  model%cts(2)%et = ts(:,8)   ! e2t

  !  Calculating the predicted values:
  !  - for mu we need g12(yt) and rt
  !  - for nu we need eta2, g22(gnut) and rt^2
  call forecast_model(model, length(2), xnew1, xnew2, forecast)
  return
end subroutine predictbtsr

subroutine gradient(n, npar, nd, D, T, h, grad)
  !----------------------------------------------------------------------------------------------
  !
  !  Subroutines to calculate the gradient matrix G
  !
  !----------------------------------------------------------------------------------------------
  !  Calculates the gradient matrix G with blocks
  !      Gr = T1 * h1 * Dr + T2 * h2 * Dr
  !      Gl = T2 * h2 * Dl
  !
  ! Input
  !   n: sample size
  !   nr: number of unknown parameters in part 1
  !   nl: number of unknown parameters in part 2
  !   nl: number of columns in D = [Dro Mrho Dlambda]
  !   T: matrix with the vectors of derivatives dmu/deta1 and dnu/deta2
  !   h: matrix with the vector of derivatives dl / mu and dl / nu
  !
  ! Output
  !   grad: the matrix with entries dl_i / dgamma_j
  !----------------------------------------------------------------------------------------------
  ! February, 2025
  !  - Added this function
  use base, only : dp
  implicit none
  integer :: n, npar(2), nd
  real(dp) :: T(n, 2)
  real(dp) :: h(n, 2)
  real(dp) :: D(n, nd)
  real(dp) :: grad(n, max(1, sum(npar)))
  integer :: j

  grad = 0.d0
  if (npar(1) > 0) then
     do j = 1, npar(1)
        ! T1 * h1 * Dr
        grad(:, j) = T(:, 1) * h(:, 1) * D(:, j)
        ! T2 * h2 * Mr
        if (npar(2) > 0) grad(:, j) = grad(:, j) + T(:, 2) * h(:, 2) * D(:, npar(1) + j)
     end do
  end if
  if (npar(2) > 0) then
     do j = 1, npar(2)
        ! T2 * h2 * Dl
        grad(:, npar(1) + j) = T(:, 2) * h(:, 2) * D(:, 2 * npar(1) + j)
     end do
  end if
  return
end subroutine gradient
