module uwarfima
  use main_mod
  use base
  implicit none

contains

  subroutine loglik_uwarfima(model, npar, par, sll, U)
    !------------------------------------------------------------------
    !
    !   Log-likelihood: UWARFIMA model
    !   Also returns the score vector if required    
    !
    !------------------------------------------------------------------
    implicit none
    integer, intent(in) :: npar
    real(dp), intent(in) :: par(npar)
    type(argsModel), intent(inout) :: model
    real(dp), intent(out) :: sll, U(npar)

    call loglik_generic(llk_uw, dllk_uw, model, npar, par, sll, U)
    if(sll < -Huge(1.d0)) sll = -Huge(1.d0)
    if(sll > Huge(1.d0)) sll = Huge(1.d0)
    return
  end subroutine loglik_uwarfima

  subroutine loglik_uwarfima_nelder(model, npar, par, sll)
    !------------------------------------------------------------------
    !
    !   Subroutine to be used in Nelder-Mead optimization subroutine
    !
    !------------------------------------------------------------------
    implicit none
    integer, intent(in) :: npar
    real(dp), intent(in) :: par(npar)
    type(argsModel), intent(inout) :: model
    real(dp), intent(out) :: sll
    real(dp) :: par_aux(npar), U(npar)

    !===================================================
    ! Back to original scale
    !===================================================
    par_aux = par
    call  transform_par(par_aux, npar, model%bounds,.true.)

    model%sco = 0
    call loglik_uwarfima(model, npar, par_aux, sll, U)
    return
  end subroutine loglik_uwarfima_nelder

end module uwarfima

subroutine uwarfimaR(n, y, gy, ystart, nreg, xreg, xstart, mu, eta, error, escale,&
     nnew, xnew, ynew, link, npar, par, fixa, alpha, fixB, flagsb, beta, p, fixar, &
     flagsar, phi, xregar,  q, fixma, flagsma, theta, fixd, d, fixnu, pdist, inf, m, &
     llk, sll, sco, U, info, K, extra, Drho, T, E, h)
  use uwarfima
  !-----------------------------------------------------------------------------
  !
  !  Wrapper to call the subroutines from R.
  !
  !-----------------------------------------------------------------------------
  ! Input:
  !--------- Time series -------------------------------------------------------
  !     n = sample size of y 
  !     y = original time series
  !     ystart = starting value for y
  !     nreg = number of regressors
  !     xreg = matrix of regressors   
  !     xstart = vector of starting values for the regressors      
  !     escale = error scale   (0 = data, 1 = predictor)   
  !     nnew = forecasting horizon (out-of-sample)
  !     xnew = matrix with the predicted values for the  regressor
  !--------- Parameters --------------------------------------------------------
  !     link = an integer corresponding to the type of link to be used
  !     npar = number of non-fixed parameters
  !     par = parameter's vector (non-fixed values). 
  !     fixa, fixd, fixnu = 0 (non-fixed) or 1 (fixed)
  !     alpha = value of alpha (if fixed). Dummy value when alpha is non-fixed.
  !     fixB, fixar, fixma = number of lags to be fixed 
  !     flagsb, flagsar, flagsma = fixed lags
  !     beta, phi, theta = fixed values for beta, phi and theta. These will be 
  !                        dummy variables if the values are non-fixed        
  !     p, q = degre of ar and ma polinomials
  !     d = value of d (if fixed). Dummy variable in case d is non-fixed.
  !     pdist = value of lambda (if fixed) and rho. Dummy value when lambda is non-fixed.
  !     inf = truncation point for infinite sums (only if d > 0)
  !--------- log-likelihood ----------------------------------------------------
  !     m = starting point for loglikelihood
  !     llk = 0 or 1 (skip or calculate the log-likelihood)
  !     sll = log-likelihood
  !     sco = 0 or 1 (skip or calculate the score)
  !     info = 0 or 1 (skip or calculate the information matrix)
  !     extra = 0 or 1 (skip or return extra matrices)
  !
  !-----------------------------------------------------------------------------      
  ! Output:      
  !-----------------------------------------------------------------------------      
  !     gy = transformed time series g2(y)   
  !     mu = conditional mean
  !     eta = g1(mu)
  !     error = y - mu or g(y) - g(mu)
  !     ynew = out-of-sample forecast
  !     U = score vector
  !     K = information matrix
  !     Drho, T, E, h = matrices used to calculate U and K
  !
  !-----------------------------------------------------------------------------
  implicit none
  ! time series
  integer :: n, nreg, link(2), inf, m, nnew, escale
  real(dp) :: y(n), gy(n), xreg(n, max(1, nreg)), xnew(max(1,nnew), max(1, nreg))
  real(dp) :: mu(n), eta(n), error(n), ynew(max(1,nnew)), ystart, xstart(max(1,nreg))
  ! parameters
  integer :: npar, fixa, fixB, fixar, fixma, fixd, fixnu, p, q , xregar  
  integer :: flagsb(max(1,fixB)), flagsar(max(1,fixar)), flagsma(max(1,fixma))
  real(dp) :: par(npar), alpha, d, pdist(2), beta(max(1,fixB))
  real(dp) :: phi(max(1,fixar)), theta(max(1,fixma))
  ! likelihood
  integer :: llk, sco, info, extra
  real(dp) :: sll, U(max(1, npar*sco))
  real(dp) :: K(max(1, npar*info), max(1, npar*info))
  real(dp) :: T(max(1, n*extra)), h(max(1, n*extra))
  real(dp) :: Drho(max(1,n*extra), max(1, (npar - 1 + fixnu)*extra))
  real(dp) :: E(max(1,n*extra), 1+2*(1-fixnu)*extra)
  !
  type(argsModel) :: model

  !---------------------------------------------------------------------------------
  ! allocating matrices and vectors and setting variables fixed values
  !--------------------------------------------------------------------------------   
  call get_model(model, n, y, 0.d0, 1.d0, ystart, gy, nreg, xreg, xstart, link,  &
       escale,0, npar, par, fixa, alpha, fixB, flagsb, beta, p, fixar, flagsar, phi, xregar,&
       q, fixma, flagsma, theta, fixd, d, fixnu, pdist(1), inf, llk, sco, info, m)
  !--------------------------------------------------       
  ! passing the value of rho for the argD argument. 
  !--------------------------------------------------
  model%argsD%arg1 = pdist(2)

  call final_model(model, npar, par, fixnu, n, mu, eta, error, nnew, nreg, &
       xnew, ynew, inf, sll, sco, U, info, K, extra, Drho, T, E, h, &
       llk_uw, dllk_uw, Ed2llk_uw)

  return
end subroutine uwarfimaR

subroutine optimnelderuwarfimaR(npar, par, nbd, lower, upper, n, y, gy, ystart, nreg, &
     xreg, xstart, mu, eta, error, escale, nnew, xnew, ynew, link, fixa, alpha, &
     fixB, flagsb, beta, p, fixar, flagsar, phi, xregar, q, fixma, flagsma, theta, fixd, d, &
     fixnu, pdist, inf, m, sll, sco, U, info, K, extra, Drho, T, E, h, iprint, stopcr, &
     maxit, neval, conv)
  use uwarfima
  implicit none
  !--------------------------------------------
  ! iprint < 0 = no print
  ! stopcr = stopping critereon  (1.d-4)
  !
  ! conv = 0 FOR SUCCESSFUL TERMINATION
  !      = 1 IF MAXIMUM NO. OF FUNCTION EVALUATIONS EXCEEDED
  !      = 2 IF NOP < 1 .or. STOPCR <= 0
  !---------------------------------------------
  integer :: n, nreg, link(2), inf, m, nnew, escale
  real(dp) :: y(n), gy(n), xreg(n, max(1, nreg)), xnew(max(1,nnew), max(1,nreg))
  real(dp) :: mu(n), eta(n), error(n), ynew(max(1,nnew)), ystart, xstart(max(1,nreg))
  integer :: npar, fixa, fixB, fixar, fixma, fixd, fixnu, p, q, xregar   
  integer :: flagsb(max(1,fixB)), flagsar(max(1,fixar)), flagsma(max(1,fixma))
  integer ::  nbd(npar), neval
  real(dp) :: par(npar), alpha, d, pdist(2), beta(max(1,fixB))
  real(dp) :: phi(max(1,fixar)), theta(max(1,fixma))
  integer ::  sco, info, extra, maxit, iprint, conv
  real(dp) :: sll, U(max(1, npar*sco)), T(max(1,n*extra)), h(max(1,n*extra))
  real(dp) :: K(max(1, npar*info),max(1, npar*info))
  real(dp) :: Drho(max(1,n*extra), max(1,npar - 1 + fixnu))
  real(dp) :: E(max(1,n*extra), 1+2*(1-fixnu)*extra)
  real(dp) :: stopcr, lower(npar), upper(npar)
  type(argsModel) :: model

  !-------------------------------------------------------------------------
  ! allocating matrices and vectors and setting variables fixed values
  ! score vector related matrices will not be allocated at this point
  !-------------------------------------------------------------------------
  call get_model(model, n, y, 0.d0, 1.d0, ystart, gy, nreg, xreg, xstart, link, &
       escale, 0, npar, par, fixa, alpha, fixB, flagsb, beta, p, fixar, flagsar, phi, xregar,&
       q, fixma, flagsma, theta, fixd, d, fixnu, pdist(1), inf, 1, 0, 0, m)
  model%argsD%arg1 = pdist(2)

  !-------------------------------------------------------------------------------
  ! Nelder-Mead
  !-------------------------------------------------------------------------------
  call optim(loglik_uwarfima_nelder, model, npar, par, nbd, lower, upper, sll, iprint,&
       stopcr, maxit, neval, conv)

  !----------------------------------------------------------------------------------
  ! Summary
  !----------------------------------------------------------------------------------
  if(sco + info > 0) call allocate_SI(model, model%SI)
  call final_model(model, npar, par, fixnu, n, mu, eta, error, nnew, nreg, &
       xnew, ynew, inf, sll, sco, U, info, K, extra, Drho, T, E, h, &
       llk_uw, dllk_uw, Ed2llk_uw)

  return      
end subroutine optimnelderuwarfimaR

subroutine optimlbfgsbuwarfimaR(npar, par, nbd, lower, upper, n, y, gy, ystart, nreg, &
     xreg, xstart, mu, eta, error, escale, nnew, xnew, ynew, link, &
     fixa, alpha, fixB, flagsb, beta, p, fixar, flagsar, phi, xregar, q, fixma, flagsma, &
     theta, fixd, d, fixnu, pdist, inf, m, sll, U, info, K, extra, Drho, T, E, h, iprint,&
     factr, pgtol, maxit, neval, conv) 
  use uwarfima
  implicit none
  integer :: n, nreg, link(2), inf, m, conv, nnew, escale
  real(dp) :: y(n), gy(n), xreg(n, max(1, nreg)), xnew(max(1,nnew), max(1,nreg))
  real(dp) :: mu(n), eta(n), error(n), ynew(max(1,nnew)), ystart, xstart(max(1,nreg))
  integer :: npar, fixa, fixB, fixar, fixma, fixd, fixnu, p, q, xregar  
  integer :: flagsb(max(1,fixB)), flagsar(max(1,fixar)), flagsma(max(1,fixma))
  integer ::  nbd(npar), neval
  real(dp) :: par(npar), alpha, d, pdist(2), beta(max(1,fixB))
  real(dp) :: phi(max(1,fixar)), theta(max(1,fixma))
  integer ::  info, extra, maxit, iprint
  real(dp) :: sll, U(npar), T(max(1,n*extra)), h(max(1,n*extra))
  real(dp) :: K(max(1, npar*info),max(1, npar*info))
  real(dp) :: Drho(max(1,n*extra), max(1,npar - 1 + fixnu))
  real(dp) :: E(max(1,n*extra), 1+2*(1-fixnu)*extra)
  real(dp) :: factr, pgtol, lower(npar), upper(npar)
  character(len = 60) :: conve
  type(argsModel) :: model
  !---------------------------------------------------------------------------------
  ! allocating matrices and vectors and setting variables fixed values
  !---------------------------------------------------------------------------------
  call get_model(model, n, y, 0.d0, 1.d0, ystart, gy, nreg, xreg, xstart, link, &
       escale, 0, npar, par, fixa, alpha, fixB, flagsb, beta, p, fixar, flagsar, phi, xregar,&
       q, fixma, flagsma, theta, fixd, d, fixnu, pdist(1), inf, 1, 1, 0, m)
  model%argsD%arg1 = pdist(2)


  !----------------------------------------------------------------------------------
  ! L-BFGS-B
  !----------------------------------------------------------------------------------
  call optim(loglik_uwarfima, model, npar, par, nbd, lower, upper, sll, U, iprint, &
       factr, pgtol, maxit, neval, conv, conve)

  !----------------------------------------------------------------------------------
  ! Summary
  !----------------------------------------------------------------------------------
  call final_model(model, npar, par, fixnu, n, mu, eta, error, nnew, nreg, &
       xnew, ynew, inf, sll, 1, U, info, K, extra, Drho, T, E, h, &
       llk_uw, dllk_uw, Ed2llk_uw)

  return
end subroutine optimlbfgsbuwarfimaR

subroutine simuwarfimaR(n, burn, pdist, alpha, nreg, beta, p, phi, &
     q, theta, d, linkg, xreg, xregar, yt, ystart, xstart, mu, eta, error, &
     escale, ns, seed, rngtype, inf, rev)
  use uwarfima
  !#---------------------------------------------------------------------
  !#
  !#  Simulating a UWARFIMA model
  !#  pdist = (lambda, rho) 
  !#  mu_t is the time varying parameter. 
  !#
  !#---------------------------------------------------------------------
  implicit none
  integer :: n, burn, nreg, p, q, linkg(2), ns, inf
  integer :: seed(ns), rngtype, rev, xregar, escale
  real(dp) :: alpha, beta(max(1,nreg)), pdist(2), d
  real(dp) :: xreg(n+burn,max(1, nreg))
  real(dp) :: ystart, xstart(max(1,nreg))
  real(dp) :: phi(max(1,p)), theta(max(1,q))
  real(dp) :: yt(n+burn), mu(n+burn), eta(n+burn), error(n+burn)
  type(argsLink) :: argsL(2)

  call sim_model(ruw, n, burn, 2, pdist, alpha, nreg, beta, p, phi, q, theta, &
       d, linkg, xreg, xregar, yt, 0.d0, 1.d0, ystart, xstart, mu, eta, &
       error, escale, ns, seed, rngtype, inf, rev, argsL)

  return
end subroutine simuwarfimaR

