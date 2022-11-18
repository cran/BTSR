module base
  !****************************************************************
  !
  !  This module contains the base subroutines that can be used
  !  by any model with the folowing structure:
  ! 
  !  Y ~ f(y | mu, sigma2)
  !  g1(mu) = a1 + X1*b1 + AR1 + MA1
  !  g2(sigma2) = a2 + X2*b2 + AR2 + MA2
  !
  !****************************************************************
  use main_mod    ! defines the special type variables 
  use RNG_mod     ! random number generation
  use specfun     ! special functions
  use lbfgsb      ! L-BFGS-B algorithm
  use Nelder_mead ! Nelder-Mead   algorithm
  implicit none  

  interface allocate_model
     module procedure allocate_model1
     module procedure allocate_model2
  end interface allocate_model

  interface start_par
     module procedure start_par1
     module procedure start_par2
  end interface start_par

  interface calc_Us
     module procedure calc_us1
     module procedure calc_us2
  end interface calc_Us

  interface calc_K
     module procedure calc_K1
     module procedure calc_K2
  end interface calc_K

  interface optim
     module procedure optim_lbfgsb
     module procedure optim_nelder
  end interface optim


contains

  !**************************************************************************************************
  !
  !                  Subroutines related to the link function
  !
  !**************************************************************************************************
  function linkfun(x, args) result(lk)
    !-------------------------------------
    !
    ! link function
    !
    !-------------------------------------
    ! outros links podem ser implementados posteriormente
    ! link, invlik e difflink devem ser coerentes
    implicit none
    type(argslink), intent(in) :: args
    real(dp), intent(in) :: x
    real(dp) :: lk, a, b

    a = args%lower
    b = args%upper
    lk = 0.d0

    select case(args%link)
    case(0)
       lk  = args%a*x
    case(1)
       lk = log((x - a)/(b - x))
    case(2)
       lk = log(x - a)
    case(3)
       lk = log(-log((x - a)/(b - a)))
    case(4)
       lk = log(-log(1.d0 - (x - a)/(b - a)))
    end select
    return   
  end function linkfun

  function linkinv(x,args) result(inv)
    !-------------------------------------
    !
    ! Inverse link function
    !
    !-------------------------------------
    implicit none
    type(argslink), intent(in) :: args
    real(dp), intent(in) :: x
    real(dp) :: inv, a, b

    a = args%lower
    b = args%upper
    inv = 0.d0

    select case(args%link)
    case(0)
       inv = x/args%a
    case(1)
       inv = b/(1.d0+exp(-x))
       if(a /= 0.d0) inv = inv + a/(exp(x) + 1.d0)
    case(2)
       inv = a + exp(x)
    case(3)
       inv =  a + (b-a)*exp(-exp(x))
    case(4)
       inv = b - (b-a)*exp(-exp(x))
    end select
    return   
  end function linkinv

  function diflink(x, args) result(dl)
    !-------------------------------------
    !
    ! link function derivative
    !
    !-------------------------------------
    implicit none
    type(argslink), intent(in) :: args    
    real(dp) :: x
    real(dp) :: dl, a, b

    a = args%lower
    b = args%upper
    dl = 0.d0

    select case(args%link)
    case(0)
       dl = args%a
    case(1)
       dl = (b-a)/((a+b)*x - x**2 - a*b)
    case(2)
       dl = 1.d0/(x-a)
    case(3)
       dl = 1.d0/((x - a)*log((x - a)/(b - a)))
    case(4)
       dl = 1.d0/((x - b)*log((b - x)/(b - a)))
    end select
    return
  end function diflink

  !**********************************************************************************************
  !
  !               Subroutines related to polinomials 
  !
  !***********************************************************************************************
  subroutine pi_f(d, inf, pik)
    !-------------------------------------------------------------
    !
    !  coeficients of (1-z)^{-d}
    !
    !-------------------------------------------------------------
    implicit none
    integer, intent(in) :: inf
    real(dp), intent(in) :: d
    real(dp), intent(out) :: pik(0:inf)
    integer :: j 
    pik = 0.d0   
    pik(0) = 1.d0
    if(d == 0.0d0) return
    do j = 1,inf
       pik(j) = pik(j-1)*(j-1+d)/j
    end do
    return
  end subroutine pi_f

  subroutine vc_f(d, theta, q, inf, ck)
    !--------------------------------------------------------
    !
    !  coeficients of theta(z)*(1-z)^{-d}
    ! 
    !  q = degree of theta
    !  inf = truncation point
    !  theta = (1, theta_1, ..., theta_q)
    !  ck(j) = sum_{i=0}^{min{j,q}} theta(i)*pik(j-i)
    !
    !-------------------------------------------------------
    implicit none
    integer, intent(in) :: inf, q
    real (dp), intent(in) :: d
    real (dp), intent(in) :: theta(0:q)
    real (dp), intent(out) :: ck(0:inf)
    real (dp) :: pik(0:inf)
    integer:: i,k

    ck = 0.d0      
    if(d == 0.0d0) then
       ! to avoid extra calculation
       ck(0:q) = theta
       return
    end if

    call pi_f(d, inf, pik)
    do k = 0,q
       do i = 0,k
          ck(k) = ck(k) + pik(k-i)*theta(i)
       end do
    end do
    do k = (q + 1),inf
       do i = 0,q
          ck(k) = ck(k) + pik(k-i)*theta(i)
       end do
    end do
    return
  end subroutine vc_f

  !*************************************************************************************************
  !
  !                Subroutines used to calculate the conditional time series
  !                associated to the model   
  !
  !**************************************************************************************************
  subroutine mu_calc(n,yt,gy,ystart,nreg,xreg,xstart,ut,eta,error,escale,&
       alpha,beta,p,phi,xregar,inf,vc,m,argsL)
    !-----------------------------------------------------------------------------------
    !
    ! Helper: Recursion for the conditional mean: ARMA-type structures
    !
    ! each particular model will have extra allocation settings if necessary
    !
    ! Input:
    !   n: sample size
    !
    !   yt: observed time series
    !
    !   gy: g^(2)(yt). Transformed time series to be used in the AR part of the model
    !
    !   ystart: starting value for yt (will raplace yt, for t < 1)
    !
    !   nreg: number of regressors
    !
    !   xreg: matrix of regressors
    !
    !   xstart: starting value for xreg (will raplace xreg, for t < 1)
    !
    !   escale: controls the scale of the error.
    !           0 = data scale. 1 = predictive scale
    !
    !   alpha, beta, phi: parameters of the model
    !
    !   p: order of AR polinomial
    !
    !   xregar: 0 = xreg is included only in the intercept
    !           1 = xreg is also included in the AR part.
    !
    !   inf: trunction point for infinite sums
    !
    !   vc: coeficients ck
    !
    !   m: starting point for the sum of likelihood
    !
    !   argsL: arguments for link functions g^(1) and g^(2)
    !
    ! Output
    !   ut: conditional time series mut
    !
    !   eta: conditional time series g^(1)(mut)
    !
    !   error: the error term rt. Depends on escale.
    !
    !------------------------------------------------------------------------------------
    implicit none
    type(argslink), intent(in) :: argsL(2)
    integer, intent(in) :: m, n, nreg, p, inf, xregar, escale
    real(dp), intent(out) :: ut(n), eta(n), error(n)
    real(dp), intent(in) :: gy(n), xreg(n,max(1,nreg)), ystart, xstart(max(1, nreg))   
    real(dp), intent(in) :: phi(max(1,p)), beta(max(1,nreg)), yt(n)
    real(dp), intent(in) :: alpha, vc(0:inf)
    integer :: t, j   
    real(dp) :: a, b
    real(dp) :: xb, ytemp, gytemp     

    ! To avoid the bounds:
    a = argsL(1)%lower
    b = argsL(1)%upper

    ! initialization 
    error = 0.d0
    eta = 0.d0
    ytemp = 0.d0
    XB = 0.d0

    ! starting values 
    if(p > 0) then
       if(a < ystart .and. ystart < b) ytemp = linkfun(ystart, argsL(2))
       if(xregar == 1 .and. nreg > 0) xb = sum(xstart*beta)
    end if

    do t = (m+1), n
       !------------------------------------
       ! eta(t) = alpha + x*b + AR + MA
       !------------------------------------
       eta(t) = alpha 
       if(nreg > 0) eta(t) = eta(t) + sum(xreg(t,:)*beta)

       do j = 1,p                    
          if(t-j > 0) then
             ! update g2(y)
             ytemp = gy(t-j)
             ! check if xreg is used in AR recursion
             if(xregar == 1 .and. nreg > 0) xb = sum(xreg(t-j,:)*beta)
          end if
          eta(t) = eta(t) + (ytemp - xb)*phi(j)          
       end do

       do j = 1, min(t-1, inf)
          ! sum(c(j)*r(t-j)), only exists if d > 0 or q > 0
          ! if inf < 1 does not enters the loop
          eta(t) =  eta(t) + vc(j)*error(t-j) 
       end do

       ut(t) = linkinv(eta(t), argsL(1))

       ! to avoid the bounds:
       if(ut(t) <= a) then
          ut(t) = a + epsilon(1.d0)
          eta(t) = linkfun(ut(t), argsL(1))
       elseif(ut(t) >= b) then
          ut(t) = b - epsilon(1.d0)
          eta(t) = linkfun(ut(t), argsL(1))
       end if

       ! calculating rt
       if(escale == 0) then
          error(t) = yt(t) - ut(t)
       else
          if(argsL(1)%link /= argsL(2)%link) then
             gytemp = linkfun(yt(t), argsL(1))
          else
             gytemp = gy(t)
          end if
          error(t) = gytemp - eta(t)
       end if
    end do

    return
  end subroutine mu_calc

  subroutine mu_forecast(model, vc, nnew, xhat, yhat)
    implicit none
    type(argsModel), intent(inout) :: model
    integer, intent(in) :: nnew
    real(dp), intent(in) :: xhat(nnew, max(1, model%beta(1)%length))
    real(dp), intent(in) :: vc(0:model%inf(1))
    real(dp), intent(out) :: yhat(nnew)
    real(dp) :: XB((model%n+1-model%ar(1)%length):(model%n+nnew))
    real(dp) :: gy((model%n+1-model%ar(1)%length):(model%n+nnew))
    real(dp) :: gyhat(nnew)
    integer :: t, i

    ! xreg*beta, t = n+1-p,...,n+1,...,n+nnew
    XB = 0.d0
    if(model%beta(1)%length > 0) then
       do t = (model%n+1-model%ar(1)%length),model%n
          XB(t) = sum(model%cts(1)%xreg(t,:)*model%beta(1)%par)
       end do
       do t = 1, nnew
          XB(model%n+t) = sum(xhat(t,:)*model%beta(1)%par)
       end do
    end if

    ! g2(y) 
    if(model%ar(1)%length > 0) gy((model%n+1-model%ar(1)%length):model%n) = &
         model%gy((model%n+1-model%ar(1)%length):model%n)  

    do t = 1,nnew
       ! ghat(y) = a + x*b
       gyhat(t) = model%alpha(1)%par(1) + XB(model%n+t)
       ! ghat(y) = a + x*b + MA
       do i = t, min(model%n+t-1, model%inf(1))  
          gyhat(t) = gyhat(t) + vc(i)*model%error(model%n+t-i)
       end do
       ! ghat(y) = a + x*b + MA + AR
       do i = 1,model%ar(1)%length   ! if p < 1 ignores the loop
          gyhat(t) = gyhat(t) + model%ar(1)%par(i)*gy(model%n+t-i)
          if(model%cts(1)%xregar == 1) gyhat(t) = gyhat(t) - model%ar(1)%par(i)*XB(model%n+t-i)
       end do
       ! yhat = g1^{-1}(gyhat)
       yhat(t) = linkinv(gyhat(t),model%argsL(1))
       ! g2(y) to be used in the AR recursion
       if(model%argsL(1)%link == model%argsL(2)%link) then
          gy(model%n+t) = gyhat(t)
       else
          gy(model%n+t) = linkfun(yhat(t),model%argsL(2))
       end if
    end do
    return
  end subroutine mu_forecast

  subroutine sigma_calc(n, error, r20, nreg, xreg, xstart, ut, ustart, eta, &
       alpha, beta, p, phi, xregar, inf, vc, m, argsL)
    !-----------------------------------------------------------------------------------
    !
    ! Helper: Recursion for the conditional mean: ARMA-type structures
    !
    ! each particular model will have extra allocation settings if necessary
    !
    ! Input:
    !   n: sample size
    !
    !   error: the error term rt corresponding to part 1 of the model
    !
    !   r20: starting values for the square of rt
    !
    !   nreg: number of regressors
    !
    !   xreg: matrix of regressors 
    !
    !   xstart: starting value for xreg (will raplace xreg, for t < 1)
    !
    !   ustart: starting values for sigmat^2 (ut)
    !
    !   alpha, beta, phi: parameters of the model
    !
    !   p: order of AR polinomial
    !
    !   xregar: 0 = xreg is included only in the intercept
    !           1 = xreg is also included in the AR part.
    !
    !   inf: trunction point for infinite sums
    !
    !   vc: coeficients ck
    !
    !   m: starting point for the sum of likelihood
    !
    !   argsL: arguments for link functions g^(1) and g^(2)
    !
    ! Output
    !   ut: conditional time series sigma_t^2
    !
    !   eta: conditional time series g^(1)(sigma_t^2)
    !
    !------------------------------------------------------------------------------------ 
    implicit none
    type(argslink), intent(in) :: argsL(2)
    integer, intent(in) :: m, n, nreg, p, inf, xregar
    real(dp), intent(in) :: error(n), r20, ustart
    real(dp), intent(in) :: xreg(n, max(1,nreg)), xstart(max(1,nreg))
    real(dp), intent(out) :: ut(n), eta(n)
    real(dp), intent(in) :: phi(max(1,p)), beta(max(1, nreg))
    real(dp), intent(in) :: alpha, vc(0:inf)    

    real(dp) :: error2(-inf:n), gs(-p:n)
    integer :: t, j   
    real(dp) :: xb, a, b

    ! To avoid the bounds:
    a = argsL(1)%lower
    b = argsL(1)%upper

    ! initializing
    error2(-inf:0) = r20   ! t < 1
    error2(1:n) = error**2 ! t > 0
    XB = 0.d0 
    eta = 0.d0  ! g^1(sigmat)
    gs = 0.d0   ! g^2(sigmat)

    if(p > 0) then
       if(abs(ustart) /= 0.d0) gs(-p:0) = linkfun(ustart, argsL(2))
       if(xregar == 1 .and. nreg > 0) xb = sum(xstart*beta)
    end if

    do t = (m+1), n
       !------------------------------------
       ! eta(t) = alpha + x*b + AR + MA
       !------------------------------------
       eta(t) = alpha 
       if(nreg > 0) eta(t) = eta(t) + sum(xreg(t,:)*beta)

       do j = 1,p                    
          if(t-j > 0) then
             ! check if xreg is used in AR recursion
             if(xregar == 1 .and. nreg > 0) xb = sum(xreg(t-j,:)*beta)
          end if
          eta(t) = eta(t) + (gs(t-j) - xb)*phi(j)          
       end do

       do j = 1, inf
          ! sum(c(j)*r2(t-j)), only exists if d > 0 or q > 0
          ! if inf < 1 does not enters the loop
          eta(t) =  eta(t) + vc(j)*error2(t-j) 
       end do

       ut(t) = linkinv(eta(t), argsL(1))

       ! to avoid the bounds:
       if(ut(t) <= a) then
          ut(t) = a + epsilon(1.d0)
          eta(t) = linkfun(ut(t), argsL(1))
       elseif(ut(t) >= b) then
          ut(t) = b - epsilon(1.d0)
          eta(t) = linkfun(ut(t), argsL(1))
       end if

       if(argsL(1)%link /= argsL(2)%link) then
          gs(t) = linkfun(ut(t), argsL(2))
       else
          gs(t) = eta(t)
       end if
    end do
    return
  end subroutine sigma_calc


  !***************************************************************************************************
  !
  !           Subroutines to calculate the generic log-likelihood
  !
  !****************************************************************************************************  
  subroutine loglik_generic(llk_dist, dllk_dist, model, npar, par, sll, U)
    !------------------------------------------------------------------
    !
    !   Log-likelihood: Generic model
    !   Also returns the score vector if required    
    !
    !------------------------------------------------------------------
    implicit none
    integer, intent(in) :: npar
    real(dp), intent(in) :: par(npar)
    type(argsModel), intent(inout) :: model
    real(dp), intent(out) :: sll, U(npar)
    real(dp) :: vc(0:model%inf(1))
    interface
       function llk_dist(m, n, y, mu, nu, argsD) result(sll)
         import :: dp, argsDist
         implicit none
         integer, intent(in) :: m, n
         real(dp), intent(in) :: y(n), mu(n), nu
         type(argsDist), intent(in) :: argsD
         real(dp) :: sll
       end function llk_dist
    end interface
    interface
       subroutine dllk_dist(m,n,y,n1,mut,skipmu,n2,nut,skipnu,dllmu,dllnu,argsD)
         import :: dp, argsDist
         implicit none
         type(argsDist), intent(in) :: argsD
         integer, intent(in) :: m, n, n1, n2, skipmu, skipnu
         real(dp), intent(in) :: y(n), mut(n1), nut(n2)
         real(dp), intent(out) :: dllmu(n*(1-skipmu)+skipmu)
         real(dp), intent(out) :: dllnu(n*(1-skipnu)+skipnu)
       end subroutine dllk_dist
    end interface

    ! Initializing parameters 
    call start_par(par, model, vc, 1)

    ! Calculating recursions that do not depend on the distribution
    call mu_calc(model%n, model%y, model%gy, model%ystart, &
         model%cts(1)%nreg,model%cts(1)%xreg,model%cts(1)%xstart, &
         model%cts(1)%ut,model%cts(1)%eta,model%error, model%escale, &
         model%alpha(1)%par(1),model%beta(1)%par,model%ar(1)%length, model%ar(1)%par,&
         model%cts(1)%xregar, model%inf(1),vc, model%m, model%argsL(1:2))

    !-----------------------------------------
    ! log-likelihood for Generic model
    !-----------------------------------------
    sll = llk_dist(model%m, model%n, model%y, model%cts(1)%ut,&
         model%nu%par(1), model%argsD)
    sll = -sll

    U = 0.d0
    if(model%sco == 0) return

    !-----------------------------------------
    !    Score vector
    !-----------------------------------------    
    call U_generic(dllk_dist, model, vc, U)
    U = -U
    return
  end subroutine loglik_generic


  !***************************************************************************************************
  !
  !           Subroutines to calculate extra matrices used to calcualte U an K
  !
  !****************************************************************************************************  
  subroutine calc_T(argsL, m, n, ut, T)
    !-----------------------------------------------------------------------------------
    !
    ! Helper: Score vector. Calculates T1
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !------------------------------------------------------------------------------------
    implicit none
    type(argsLink), intent(inout) :: argsL
    integer, intent(in) ::  n, m
    real(dp), intent(in) ::  ut(n)
    real(dp), intent(out) :: T(n) 
    integer :: i
    T = 0.d0
    do i = (m +1), n
       T(i) = 1.d0/diflink(ut(i), argsL)
    end do
    return
  end subroutine calc_T


  !********************************************************************************************
  !
  !   Subroutines to calculate  and organize the derivatives deta_t/dgamma_j
  !   and the associated matrices Drho and Dlambda
  !
  !*******************************************************************************************  
  function digamma(x) result(fn)
    !------------------------------------
    ! digamma function - psi(x)
    !
    ! Implemented by Taiane S. Prass
    ! March, 2018
    ! Calculates the value of psi(x)
    !-------------------------------------
    implicit none
    real(dp), intent(in) :: x
    real (dp) :: fn
    fn = psi(x)
  end function digamma

  subroutine deta1_drho(model, SI, vc)
    !-----------------------------------------------------------------------------------
    !
    ! Helper: Recursion for the score vector. Calculates deta1_drho
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !------------------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(in) :: model
    type(argsSI), intent(inout) :: SI
    real(dp), intent(in) :: vc(0:model%inf(1))
    real(dp) :: diff(model%n)
    real(dp), allocatable :: diffm(:,:)
    integer :: i, j, jj,  k, t, n1, n2, es   
    real(dp) :: soma, pik(0:model%inf(1)), ies0

    call pi_f(model%d(1)%par(1), model%inf(1), pik)   

    es = model%escale
    ies0 = 1.d0

    ! nu: will be calculated outside this subroutine

    n1 = 1
    n2 = model%alpha(1)%fit    
    ! alpha
    if(n2 >= n1) then       
       !
       ! diff(t) = 1 - sum(ck*diff(t-k))
       !
       diff = 0.d0
       do t = (model%m +1),model%n
          diff(t) = 1.d0
          do j = 1, min(t-1, model%inf(1))
             if(es == 0) ies0 = SI%T(1)%par(t-j)
             diff(t) = diff(t) - vc(j)*ies0*diff(t-j)
          end do
       end do
       SI%deta(1,1)%dalpha(:,1) = diff
    end if

    n1 = n2 + 1
    n2 = n2 + model%beta(1)%fit
    ! betas    
    if(n2 >= n1) then
       !
       ! diff(t) = x(t,j) - sum(phi(k)*x(t-k,j)) - sum(c(k)*diff(t-k))
       !
       call safe_allocate(diffm, model%n,model%beta(1)%fit)
       diffm = 0.d0
       do j = 1, model%beta(1)%fit
          jj = model%beta(1)%lags(j)  ! j-th non-fixed lag
          do t = (model%m +1), model%n
             diffm(t,j) = model%cts(1)%xreg(t,jj)
             ! if p = 0, then ar = 0
             if(model%cts(1)%xregar == 1) then
                do k = 1, min(t-1, model%ar(1)%length)
                   diffm(t,j) = diffm(t,j) - model%ar(1)%par(k)*model%cts(1)%xreg(t-k,jj)
                end do
             end if
             do k = 1, min(t-1, model%inf(1))
                if(es == 0) ies0 = SI%T(1)%par(t-k)
                diffm(t,j) = diffm(t,j) - vc(k)*ies0*diffm(t-k,jj)
             end do
          end do
       end do
       SI%deta(1,1)%dbeta = diffm
    end if

    n1 = n2 + 1
    n2 = n2 + model%ar(1)%fit
    ! phis
    if(n2 >= n1) then
       !
       ! diff(t) = g(y(t-j)) - X(t-j)*b -  sum(ck*diff(t-k))
       !
       call safe_allocate(diffm,model%n,model%ar(1)%fit)
       diffm = 0.d0
       do j = 1,  model%ar(1)%fit
          jj = model%ar(1)%lags(j)  ! j-th non-fixed lag
          do t =  max((jj + 1),(model%m +1)), model%n   ! so that t - jj > 0             
             diffm(t,j) =  model%gy(t-jj)
             if(model%cts(1)%xregar == 1) then
                do i = 1, model%cts(1)%nreg  ! if nreg = 0 does not enter the loop
                   diffm(t,j) = diffm(t,j) - model%cts(1)%xreg(t-jj,i)*model%beta(1)%par(i)
                end do
             end if
             do k = 1,min(t-1,model%inf(1))
                if(es == 0) ies0 = SI%T(1)%par(t-k)
                diffm(t,j) = diffm(t,j) - vc(k)*ies0*diffm(t-k,j)
             end do
          end do
       end do
       SI%deta(1,1)%dphi = diffm
    end if

    n1 = n2 + 1
    n2 = n2 + model%ma(1)%fit
    ! thetas
    if(n2 >= n1) then
       !
       ! diff(t) = sum(pik(k-s)*r(t-k), k = s,...,infty) -  sum(ck*diff(t-k))
       !
       call safe_allocate(diffm, model%n, model%ma(1)%fit)
       diffm = 0.d0
       do j = 1,model%ma(1)%fit  !if q > 0, then vc(k) is different from zero for some k =1,...,inf
          jj = model%ma(1)%lags(j)  ! j-th non-fixed lag
          do t = (model%m +1),model%n
             do k = 1, min(t-1, (jj - 1))
                if(es == 0) ies0 = SI%T(1)%par(t-k)
                diffm(t,j) = diffm(t,j) - vc(k)*ies0*diffm(t-k,j)
             end do
             do k = jj, min(t-1,model%inf(1))
                if(es == 0) ies0 = SI%T(1)%par(t-k)
                diffm(t,j) =  diffm(t,j) +  pik(k-jj)*model%error(t-k) - vc(k)*ies0*diffm(t-k,j)          
             end do
          end do
       end do
       SI%deta(1,1)%dtheta = diffm(1:model%n,:)
    end if

    n1 = n2 + 1
    n2 = n2 + model%d(1)%fit
    ! d
    if(n2 >= n1) then
       !
       ! diff(t) = sum(r(t-k)*s(k)) -  sum(ck*diff(t-k)),
       !  where s(k) = sum(theta(i)*pi(k-i)*digamas(i))
       !
       diff = 0.d0
       do t = (model%m +1),model%n
          do k = 1, min(t-1,model%inf(1))
             ! i = 0
             soma = pik(k)*(digamma(model%d(1)%par(1)+k)-digamma(model%d(1)%par(1)))
             do i = 1, min(k, model%ma(1)%length) !    ! theta(0) = 1. Se q = 0 nÃ£o entra no loop
                soma = soma + model%ma(1)%par(i)*pik(k-i)*(digamma(model%d(1)%par(1)+k-i) - &
                     digamma(model%d(1)%par(1)))
             end do
             if(es == 0) ies0 = SI%T(1)%par(t-k)
             diff(t) = diff(t) + model%error(t-k)*soma - vc(k)*ies0*diff(t-k)
          end do
       end do
       SI%deta(1,1)%dd(:,1) = diff
    end if

    return
  end subroutine deta1_drho

  subroutine deta2_drho(model, SI, vc)          
    !-----------------------------------------------------------------------------------
    !
    ! Helper: Recursion for the score vector. Calculates deta2/drho
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !-----------------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(in) :: model
    type(argsSI), intent(inout) :: SI
    real(dp), intent(in) :: vc(0:model%inf(2))
    real(dp) :: diff(model%n)
    real(dp), allocatable :: diffm(:,:), es, ies0
    integer :: j, k, t, n1, n2

    es = model%escale
    ies0 = 1.d0

    !
    !  diff2(t) = sum(phi(k)*diff2(t-k)) - 2*sum(c(k)*r(t-k)*diff1(t-k))
    !

    n1 = 1
    n2 = model%alpha(1)%fit    
    ! alpha
    if(n2 >= n1) then
       diff =  0.d0
       do t = 1,model%n
          do k = 1, min(t-1,model%ar(2)%length)
             diff(t) = diff(t) + model%ar(2)%par(k)*diff(t-k)
          end do
          do k = 1, min(t-1, model%inf(2))
             if(es == 0) ies0 = SI%T(1)%par(t-k)
             diff(t) = diff(t) - 2*vc(k)*model%error(t-k)*ies0*SI%deta(1,1)%dalpha(t-k,1)  
          end do
       end do
       SI%deta(2,1)%dalpha(:,1) = diff
    end if

    n1 = n2 + 1
    n2 = n2 + model%beta(1)%fit
    ! betas
    if(n2 >= n1) then
       call safe_allocate(diffm,model%n,model%beta(1)%fit)
       diffm = 0.d0
       do j = 1,model%beta(1)%fit          
          do t = 1, model%n
             do k = 1, min(t-1,model%ar(2)%length)
                diffm(t,j) = diffm(t,j) + model%ar(2)%par(k)*diffm(t-k,j)
             end do
             do k = 1, min(t-1, model%inf(2))
                if(es == 0) ies0 = SI%T(1)%par(t-k)
                diffm(t,j) = diffm(t,j) - 2*vc(k)*model%error(t-k)*ies0*SI%deta(1,1)%dbeta(t-k,j)  
             end do
          end do
       end do
       SI%deta(2,1)%dbeta = diffm
    end if

    n1 = n2 + 1
    n2 = n2 + model%ar(1)%fit
    ! phis
    if(n2 >= n1) then
       call safe_allocate(diffm, model%n,model%ar(1)%fit)
       diffm = 0.d0
       do j = 1,model%ar(1)%fit          
          do t = 1, model%n
             do k = 1, min(t-1,model%ar(2)%length)
                diffm(t,j) = diffm(t,j) + model%ar(2)%par(k)*diffm(t-k,j)
             end do
             do k = 1, min(t-1, model%inf(2))
                if(es == 0) ies0 = SI%T(1)%par(t-k)
                diffm(t,j) = diffm(t,j) - 2*vc(k)*model%error(t-k)*ies0*SI%deta(1,1)%dphi(t-k,j)  
             end do
          end do
       end do
       SI%deta(2,1)%dphi = diffm
    end if

    n1 = n2 + 1
    n2 = n2 + model%beta(1)%fit
    ! thetas
    if(n2 >= n1) then
       call safe_allocate(diffm, model%n, model%ma(1)%fit)
       diffm = 0.d0
       do j = 1,model%ma(1)%fit          
          do t = 1, model%n
             do k = 1,min(t-1,model%ar(2)%length)
                diffm(t,j) = diffm(t,j) + model%ar(2)%par(k)*diffm(t-k,j)
             end do
             do k = 1, min(t-1, model%inf(2))
                if(es == 0) ies0 = SI%T(1)%par(t-k)
                diffm(t,j) = diffm(t,j) - 2*vc(k)*model%error(t-k)*ies0*SI%deta(1,1)%dtheta(t-k,j)  
             end do
          end do
       end do
       SI%deta(2,1)%dtheta = diffm
    end if

    n1 = n2 + 1
    n2 = n2 +  model%d(1)%fit    
    ! d
    if(n2 >= n1) then
       diff =  0.d0
       do t = 1,model%n
          do k = 1, min(t-1,model%ar(2)%length)
             diff(t) = diff(t) + model%ar(2)%par(k)*diff(t-k)
          end do
          do k = 1, min(t-1, model%inf(2))
             if(es == 0) ies0 = SI%T(1)%par(t-k)
             diff(t) = diff(t) - 2*vc(k)*model%error(t-k)*SI%deta(1,1)%dd(t-k,1)  
          end do
       end do
       SI%deta(2,1)%dd(:,1) = diff
    end if

    return
  end subroutine deta2_drho

  subroutine deta2_dlambda(model, SI)
    !-----------------------------------------------------------------------------------
    !
    ! Helper: Recursion for the score vector. Calculates deta2/dlambda
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !-----------------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(in) :: model
    type(argsSI), intent(inout) :: SI
    real(dp) :: diff(model%n), pik(0:model%inf(2))
    real(dp) :: error2(-model%inf(2):model%n), soma
    real(dp), allocatable :: diffm(:,:)
    integer :: i, j, jj,  k, t, n1, n2  

    error2(-model%inf(2):0) = model%r20
    error2(1:model%n) = model%error(1:model%n)**2

    call pi_f(model%d(2)%par(1),model%inf(2), pik)   

    ! nu: will be calculated outside this subroutine

    n1 = 1
    n2 = model%alpha(2)%fit    
    ! alpha
    if(n2 >= n1) then       
       !
       ! diff(t) = 1 + sum(phi(k)*diff(t-k))
       !
       diff = 0.d0
       do t = 1,model%n
          diff(t) = 1.d0
          do k = 1, min(t-1, model%ar(2)%length)
             diff(t) = diff(t) + model%ar(2)%par(k)*diff(t-k)
          end do
       end do
       SI%deta(2,2)%dalpha(:,1) = diff
    end if

    n1 = n2 + 1
    n2 = n2 + model%beta(2)%fit
    ! betas    
    if(n2 >= n1) then
       !
       ! diff(t,j) = x(t,j) - sum(phi(k)*[diff(t-k,j) - x(t-k,j)])
       !
       call safe_allocate(diffm, model%n, model%beta(2)%fit)
       diffm = 0.d0
       do j = 1, model%beta(2)%fit
          jj = model%beta(2)%lags(j)  ! j-th non-fixed lag
          do t = 1, model%n
             diffm(t,j) = model%cts(2)%xreg(t,jj)
             ! if p = 0, then ar = 0 
             do k = 1, min(t-1, model%ar(2)%length)
                diffm(t,j) = diffm(t,j) + model%ar(2)%par(k)*diffm(t-k,j)
                if(model%cts(2)%xregar == 1) diffm(t,j) = diffm(t,j) - model%ar(2)%par(k)*model%cts(2)%xreg(t-k,jj) 
             end do
          end do
       end do
       SI%deta(2,2)%dbeta = diffm
    end if

    n1 = n2 + 1
    n2 = n2 + model%ar(2)%fit
    ! phis
    if(n2 >= n1) then
       !
       ! diff(t) = g(s2(t-j)) - X(t-j)*b +  sum(phi(k)*diff(t-k,j))
       !
       call safe_allocate(diffm, model%n, model%ar(2)%fit)
       diffm = 0.d0
       do j = 1, model%ar(2)%fit
          jj = model%ar(2)%lags(j)  ! j-th non-fixed lag
          do t =  (jj + 1), model%n   ! so that t - jj > 0             
             diffm(t,j) =  model%cts(2)%eta(t-jj)
             if(model%cts(2)%xregar == 1) then
                do i = 1, model%cts(2)%nreg
                   diffm(t,j) = diffm(t,j) - model%cts(2)%xreg(t-jj,i)*model%beta(2)%par(i)
                end do
             end if
             do k = 1, min(t-1, model%ar(2)%length)
                diffm(t,j) = diffm(t,j) + model%ar(2)%par(k)*diffm(t-k,j)
             end do
          end do
       end do
       SI%deta(2,2)%dphi = diffm         
    end if

    n1 = n2 + 1
    n2 = n2 + model%ma(2)%fit
    ! thetas
    if(n2 >= n1) then
       !
       ! diff(t,j) = r2(t-j) +  sum(phi(k)*diff(t-k,j))
       !
       call safe_allocate(diffm,-model%inf(2), model%n, 1, model%ma(2)%fit)
       diffm = 0.d0
       do j = 1,  model%ma(2)%fit
          jj = model%ma(2)%lags(j)  ! j-th non-fixed lag
          do t =  1, model%n   
             diffm(t,j) =  error2(t-jj) 
             do k = 1, min(t-1,model%ar(2)%length)
                diffm(t,j) = diffm(t,j) + model%ar(2)%par(k)*diffm(t-k,j)
             end do
          end do
       end do
       SI%deta(2,2)%dtheta = diffm            
    end if

    n1 = n2 + 1
    n2 = n2 + model%d(2)%fit
    ! d
    if(n2 >= n1) then
       !
       ! diff(t) =  sum(phi(k)*diff(t-k,j)) +  sum(r2(t-k)*s(k)),
       ! where s(k) = sum(theta(i)*pi(k-i)*digamas(i))
       !
       diff = 0.d0
       do t = 1,model%n
          do k = 1,min(t-1, model%ar(2)%length)
             diff(t) = diff(t) + model%ar(2)%par(k)*diff(t-k)
          end do
          do k = 1,model%inf(2)             
             soma = pik(k)*(digamma(model%d(2)%par(1)+k)-digamma(model%d(2)%par(1)))  
             do i = 1, model%ma(2)%length ! theta(0) = 1. Se q = 0 nÃ£o entra no loop
                soma = soma + model%ma(2)%par(i)*pik(k-i)*(digamma(model%d(2)%par(1)+k-i) - &
                     digamma(model%d(2)%par(1)))
             end do
             diff(t) = diff(t) + error2(t-k)*soma
          end do
       end do
       SI%deta(2,2)%dd(:,1) = diff
    end if

    return
  end subroutine deta2_dlambda

  subroutine fill_D(SI, fita, fitb, fitar, fitma, fitd, n, nd, D)
    !----------------------------------------------------------------
    ! Helper: creates a matrix with the derivatives deta1/drho
    !
    ! used to return values when calling from R.
    ! also used to calculate the matrix K.
    !----------------------------------------------------------------
    implicit none
    type(argsSI), intent(in) :: SI
    integer, intent(in) :: n, nd, fita, fitb, fitar, fitma, fitd
    real(dp), intent(inout) :: D(n,nd)
    integer :: n1, n2

    n1 = 1
    n2 = fita
    ! alpha
    if(n2 >= n1) D(:,n1:n2) = SI%deta(1,1)%dalpha
    n1 = n2 + 1
    n2 = n2 + fitb
    ! beta
    if(n2 >= n1) D(:,n1:n2) = SI%deta(1,1)%dbeta
    n1 = n2 + 1
    n2 = n2 + fitar
    ! phi
    if(n2 >= n1) D(:,n1:n2) = SI%deta(1,1)%dphi
    n1 = n2 + 1
    n2 = n2 + fitma
    ! theta
    if(n2 >= n1) D(:,n1:n2) = SI%deta(1,1)%dtheta
    n1 = n2 + 1
    n2 = n2 + fitd
    ! d
    if(n2 >= n1) D(:,n1:n2) = SI%deta(1,1)%dd

    return       
  end subroutine fill_D

  !**********************************************************************
  !
  !  Subroutines used to calculate product and sum of matrices 
  !  that appear in the score U and the information matrix K
  ! 
  !***********************************************************************    
  subroutine ATh(nra, nca, A, Th, P)
    !----------------------------------------------------------------------------------
    !
    ! Helper: Score vector. Calculates the product A*T(part)*h(part),
    ! where A = Drho, Dlambda or Mrho
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !------------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: nra, nca
    real(dp), intent(in) ::  A(nra, nca), Th(nra)
    real(dp), intent(out) :: P(nca)
    P = matmul(transpose(A),Th)
    return    
  end subroutine ATh

  subroutine calc_DTh(model, SI, part, Th)
    !-----------------------------------------------------------------------------------
    !
    ! Helper: Score vector. Calculates Drho and Dlambda
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !------------------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(in) :: model
    type(argsSI), intent(inout) :: SI
    integer, intent(in) :: part
    real(dp), intent(in) :: Th(model%n)

    ! alpha
    if(model%alpha(part)%fit == 1) call ATh(model%n,1,SI%deta(part,part)%dalpha,Th,SI%U(part)%Ualpha)

    ! betas    
    if(model%beta(part)%fit > 0) call ATh(model%n, model%beta(part)%fit, SI%deta(part,part)%dbeta,Th,SI%U(part)%ubeta)   

    ! phis
    if(model%ar(part)%fit > 0) call ATh(model%n, model%ar(part)%fit,SI%deta(part,part)%dphi,Th,SI%U(part)%uphi)   

    ! thetas
    if(model%ma(part)%fit > 0) call ATh(model%n, model%ma(part)%fit,SI%deta(part,part)%dtheta,Th,SI%U(part)%utheta)   

    ! d
    if(model%d(part)%fit == 1) call ATh(model%n, 1, SI%deta(part, part)%dd,Th,SI%U(part)%ud)   
    return
  end subroutine Calc_DTh

  subroutine AddM(model, SI, Th)
    !-----------------------------------------------------------------------------------
    !
    ! Helper: Score vector. Calculates Mrho (Mlambda = 0)
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !------------------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(in) :: model
    type(argsSI), intent(inout) :: SI
    real(dp), intent(in) :: Th(model%n)
    real(dp), allocatable :: dummy(:)
    real(dp) :: dm(1)

    ! alpha
    if(model%alpha(1)%fit  == 1) then
       call ATh(model%n, 1, SI%deta(2,1)%dalpha, Th, dm)
       SI%U(1)%Ualpha = SI%U(1)%Ualpha + dm
    end if

    ! betas    
    if(model%beta(1)%fit > 0) then
       call safe_allocate(dummy,model%beta(1)%fit)
       call ATh(model%n, model%beta(1)%fit, SI%deta(2,1)%dbeta, Th, dummy)   
       SI%U(1)%ubeta = SI%U(1)%ubeta + dummy         
    end if

    ! phis
    if(model%ar(1)%fit > 0) then
       call safe_allocate(dummy, model%ar(1)%fit)
       call ATh(model%n, model%ar(1)%fit, SI%deta(2,1)%dphi, Th, dummy)   
       SI%U(1)%uphi = SI%U(1)%uphi  + dummy
    end if

    ! thetas
    if(model%ma(1)%fit > 0) then
       call safe_allocate(dummy, model%ma(1)%fit)
       call ATh(model%n, model%ma(1)%fit, SI%deta(2,1)%dtheta, Th, dummy)       
       SI%U(1)%utheta = SI%U(1)%utheta + dummy
    end if

    ! d
    if(model%d(1)%fit == 1) then
       call ATh(model%n, 1, SI%deta(2,1)%dd, Th, dm)
       SI%U(1)%ud  = SI%U(1)%ud + dm
    end if

    if(allocated(dummy)) deallocate(dummy)
    return
  end subroutine AddM

  !******************************************************************
  !
  !  Subroutines used to calculate the score vectors and U_gamma and
  !  the vector of all scores U(gamma).
  !
  !  U_nu depends on the distribution!
  !
  !******************************************************************
  subroutine Unuh_dist(dllk_dist,argsD,m,n,y,mu,nu,fitnu,npar,Unu,h)
    !-------------------------------------------------------------------------
    !
    !  Generic subroutine to be called by the user
    !  Calculates U(nu) = dl/dnu and h1 = dl/dmu for regression models
    !
    !-------------------------------------------------------------------------
    implicit none
    type(argsDist), intent(in) :: argsD
    integer, intent(in) :: m, n, fitnu, npar
    real(dp), intent(inout) :: Unu, h(n)    
    real(dp), intent(in) :: y(n), mu(n), nu
    interface
       subroutine dllk_dist(m,n,y,n1,mut,skipmu,n2,nut,skipnu,dllmu,dllnu,argsD)
         import :: dp, argsDist
         implicit none
         type(argsDist), intent(in) :: argsD
         integer, intent(in) :: m, n, n1, n2, skipmu, skipnu
         real(dp), intent(in) :: y(n), mut(n1), nut(n2)
         real(dp), intent(out) :: dllmu(n*(1-skipmu)+skipmu)
         real(dp), intent(out) :: dllnu(n*(1-skipnu)+skipnu)
       end subroutine dllk_dist
    end interface

    integer :: skn, skm
    real(dp), allocatable :: dllmu(:), dllnu(:)

    !-------------------------------------------------------
    ! score for nu:  dl/dnu  (skip if nu is not estimated)
    !-------------------------------------------------------
    skn = 0
    if(fitnu == 0) skn = 1
    call safe_allocate(dllnu, max(1,n*(1-skn)))

    !-------------------------------------------------------------
    ! vector h = dl/dmu  (skip if nu is the only parameter)
    !--------------------------------------------------------------
    skm = 0
    if(npar == 1 .and. fitnu == 1) skm = 1
    call safe_allocate(dllmu, max(1,n*(1-skm)))    

    call dllk_dist(m,n,y,n,mu,skm,1,(/nu/),skn,dllmu,dllnu,argsD)

    !-------------------------------------------------------
    ! if nu is not estimated skip the calculation again  
    !-------------------------------------------------------
    if(fitnu == 0)  goto 200
    Unu = sum(dllnu)

200 continue    
    !-------------------------------------------------------------
    ! avoiding extra calculations when nu is the only parameter
    !--------------------------------------------------------------
    if(npar == 1 .and. fitnu == 1) goto 300
    ! calculates the vector h
    h = dllmu    

300 continue
    return
  end subroutine Unuh_dist

  subroutine calc_Us1(model, SI, vc)
    !-----------------------------------------------------------------------------------
    !
    ! Score vector. Calculates U(gamma)
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !------------------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    type(argsSI), intent(inout) :: SI
    real(dp), intent(in) :: vc(0:model%inf(1))
    real(dp) :: Th(model%n)
    integer :: i

    call deta1_drho(model, SI, vc)
    Th = 0.d0
    do i = (model%m +1),model%n
       Th(i) = SI%T(1)%par(i)*SI%h(1)%par(i)
    end do

    call calc_DTh(model,SI,1,Th)

    return
  end subroutine calc_Us1

  subroutine calc_Us2(model, SI, vc1, vc2)
    !-----------------------------------------------------------------------------------
    !
    ! Score vector. Calculates U(gamma)
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !------------------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    type(argsSI), intent(inout) :: SI
    real(dp), intent(in) :: vc1(0:model%inf(1))
    real(dp), intent(in) :: vc2(0:model%inf(2))
    real(dp) :: Th(model%n)
    integer :: i
    ! Drho*T1*h1
    call calc_Us1(model, SI, vc1)

    do i = 1,model%n
       Th(i) = SI%T(2)%par(i)*SI%h(2)%par(i)
    end do

    ! Drho*T1*h1 + Mrho*T2*h2
    call deta2_drho(model, SI, vc2)    
    call AddM(model, SI, Th)

    ! Dlambda*T2*h2
    call deta2_dlambda(model, SI)
    call calc_DTh(model,SI,2,Th)
    return    
  end subroutine calc_Us2

  subroutine fill_U(SI, fita, fitb, fitar, fitma, fitd, fitnu, npar, U)
    !----------------------------------------------------------------
    ! Helper: creates a vector U(gamma) with the derivatives
    ! dl/dgamma_j. For now, only implemented for the part 1 of the
    ! model.
    !
    ! used to return values when calling from R.
    !----------------------------------------------------------------
    implicit none
    type(argsSI), intent(in) :: SI
    integer, intent(in) :: npar, fita, fitb, fitar, fitma, fitd, fitnu
    real(dp), intent(inout) :: U(npar)
    integer :: n1, n2

    n1 = 1
    n2 = fita
    ! alpha
    if(n2 >= n1) U(n1:n2) = SI%U(1)%Ualpha
    n1 = n2 + 1
    n2 = n2 + fitb
    ! beta
    if(n2 >= n1) U(n1:n2) = SI%U(1)%Ubeta
    n1 = n2 + 1
    n2 = n2 + fitar
    ! phi
    if(n2 >= n1) U(n1:n2) = SI%U(1)%Uphi
    n1 = n2 + 1
    n2 = n2 + fitma
    ! theta
    if(n2 >= n1) U(n1:n2) = SI%U(1)%Utheta
    n1 = n2 + 1
    n2 = n2 + fitd
    ! d
    if(n2 >= n1) U(n1:n2) = SI%U(1)%Ud    
    n1 = n2 + 1
    n2 = n2 + fitnu
    ! nu
    if(n2 >= n1) U(n1:n2) = SI%U(1)%Unu   

    return       
  end subroutine fill_U

  subroutine U_generic(dllk_dist, model, vc, U)
    !----------------------------------------------------
    !
    !  Subroutine used to calculate the score vector
    !  for a generic model
    !
    !----------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    real(dp), intent(in) :: vc(0:model%inf(1))
    real(dp), intent(out) :: U(model%npar(1))
    interface
       subroutine dllk_dist(m,n,y,n1,mut,skipmu,n2,nut,skipnu,dllmu,dllnu,argsD)
         import :: dp, argsDist
         implicit none
         type(argsDist), intent(in) :: argsD
         integer, intent(in) :: m, n, n1, n2, skipmu, skipnu
         real(dp), intent(in) :: y(n), mut(n1), nut(n2)
         real(dp), intent(out) :: dllmu(n*(1-skipmu)+skipmu)
         real(dp), intent(out) :: dllnu(n*(1-skipnu)+skipnu)
       end subroutine dllk_dist
    end interface

    ! calculating U(nu) and h1
    call Unuh_dist(dllk_dist, model%argsD, model%m, model%n,model%y,model%cts(1)%ut, &
         model%nu%par(1), model%nu%fit, model%npar(1), model%SI%U(1)%Unu(1), model%SI%h(1)%par)

    ! calculates T1: used in U(rho)
    call calc_T(model%argsL(1), model%m, model%n, model%cts(1)%ut, model%SI%T(1)%par)

    ! score for other parameters 
    call calc_Us(model, model%SI, vc)

    ! Score vector
    call fill_U(model%SI, model%alpha(1)%fit, model%beta(1)%fit, model%ar(1)%fit, &
         model%ma(1)%fit, model%d(1)%fit, model%nu%fit, model%npar(1), U)  
    return    
  end subroutine U_generic

  !******************************************************************
  !
  !  Subroutines to calculate the information matrix K
  !
  !******************************************************************  
  subroutine calc_K1(n, T, nce, E, ncd, D, npar, K, part)
    !--------------------------------------------------------------------------------------
    !  This subroutine is called only if:
    !    - either mu or nu is not time varying. 
    !    - mu or nu is not the only parameter in the model that was estimated. 
    !
    !  Here E = (E_mu, E_(mu,nu), E_nu) for a full model.
    !  ne = 1 indicates that nu (or mu) was not fitted.
    !--------------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n, nce, ncd, npar, part
    real(dp), intent(in)  :: T(n), E(n,nce)
    real(dp), intent(in) :: D(n, ncd)
    real(dp), intent(out) :: K(npar, npar)
    integer :: i, j, il, ic

    ! K_{rho,rho} or K_{lambda,lambda}
    ic = 1
    if(part == 2) ic = nce
    do j = part, (part + ncd - 1)
       do i = part, j
          K(i,j) = sum(D(:,j)*T**2*E(:,ic)*D(:,i))
          K(j,i) = K(i,j)
       end do
    end do

    if(npar == ncd) return

    ! K_{nu, rho} or K_{mu,lambda}
    il = 1
    if(part == 1) il = ncd + 1
    do j = part, (part + ncd - 1)
       K(il,j) = sum(D(:,j)*T*E(:,2))
       K(j,il) = K(il,j)
    end do

    ! K_{nu,nu} or K_{mu,mu}    
    K(il,il) = sum(E(:,3))

    return
  end subroutine calc_K1

  subroutine calc_K2(n, nr, nl, T1, T2, E, Dr, Dl, Mr, K)
    !-----------------------------------------------------------------------------------
    !
    ! Information Matrix:  the case where mu and nu are time varying
    !
    ! each particular model will have extra allocation settings if necessary
    !
    ! E = (E_mu, E_(mu,nu), E_nu) for a full model
    !------------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n, nr, nl
    real(dp), intent(in), dimension(n)  :: T1, T2, E(n,3)
    real(dp), intent(in) :: Dr(n, nr), Dl(n, max(1,nl)), Mr(n, nr)
    real(dp), intent(out) :: K(nr+nl, nr+nl)
    integer :: j, i

    ! K_{rho,rho}
    do j = 1, nr
       do i = 1,j
          K(i,j) = sum((Dr(:,j)*T1**2*E(:,1) + Mr(:,j)*T1*T2*E(:,2))*Dr(:,i) + &
               (Dr(:,j)*T1*T2*E(:,2) + Mr(:,j)*T2**2*E(:,3))*Mr(:,i)) 
          K(j,i) = K(i,j)
       end do
    end do

    ! K_{rho,lambda}
    do j = 1, nr
       do i = (nr+1),(nr+nl)
          K(i,j) = sum((Dr(:,j)*T1*T2*E(:,2) + Mr(:,j)*T2**2*E(:,3))*Dl(:,i-nr)) 
          K(j,i) = K(i,j)
       end do
    end do

    ! K_{lambda,lambda}
    do j = (nr+1), (nr + nl)
       do i = (nr+1),j
          K(i,j) = sum(Dl(:,j-nr)*T2**2*E(:,3)*Dl(:,i-nr))
          K(j,i) = K(i,j)
       end do
    end do

    return   
  end subroutine calc_K2

  subroutine K_generic(Ed2llk_dist, SI, mu, fita, fitb, fitar, fitma, fitd, fitnu, &
       nu, m, n, npar, K, argsD)
    !
    ! Requires: mut, Drho, T
    !
    implicit none
    type(argsSI), intent(inout) :: SI
    integer, intent(in) :: fita, fitb, fitar, fitma, fitd, fitnu
    integer, intent(in) :: npar, n, m
    real(dp), intent(in) :: mu(n), nu
    real(dp), intent(out) :: K(npar, npar)
    type(argsDist), intent(in) :: argsD
    interface
       subroutine Ed2llk_dist(m,n,n1,mut,i1,n2,nut,i2,E,argsD)
         import :: dp, argsDist
         implicit none
         integer, intent(in) :: m, n, n1, n2, i1, i2
         real(dp), intent(in) :: mut(n1), nut(n2)
         real(dp), intent(out) :: E(n,max(1,i1+i2+i1*i2))
         type(argsDist), intent(in) :: argsD
       end subroutine Ed2llk_dist
    end interface
    real(dp) :: Drho(n, max(1,npar - fitnu))
    integer :: i1

    i1 = 1
    if(fitnu == npar) i1 = 0
    call safe_allocate(SI%E, n, 1+2*fitnu*i1)
    call Ed2llk_dist(m,n,n,mu,i1,1,(/nu/),fitnu,SI%E,argsD)

    if(npar == fitnu) then
       !------------------------------------------------------------
       ! nu is the only parameter to be estimated
       !------------------------------------------------------------       
       K = sum(SI%E((m+1):n,1))
       return
    end if

    !-------------
    ! Drho
    !------------
    call fill_D(SI, fita, fitb, fitar, fitma, fitd, n, npar-fitnu, Drho)
    call calc_K(n, SI%T(1)%par, 1+2*fitnu, SI%E, npar-fitnu, Drho, npar, K, 1)
    return
  end subroutine K_generic


  !*********************************************************************
  !
  !                 Optimization subroutines
  !
  !  Available methods
  !  - Nelder-Mead : requires loglik(model, npar, par)   
  !  - L-BFGS-B : requires loglik(model, npar, par, sll, score)
  !
  !********************************************************************* 
  subroutine optim_nelder(loglik, model, npar, par, nbd, lower, upper,&
       sll, iprint, stopcr, maxit, neval, conv)
    !--------------------------------------------
    ! iprint < 0 = no print
    ! stopcr = stopping critereon  (1.d-4)
    !
    ! conv = 0 FOR SUCCESSFUL TERMINATION
    !      = 1 IF MAXIMUM NO. OF FUNCTION EVALUATIONS EXCEEDED
    !      = 2 IF NOP < 1 .or. STOPCR <= 0
    !---------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    integer, intent(in) :: iprint, maxit, npar
    integer, intent(inout) :: nbd(npar)
    integer, intent(out) ::  neval, conv
    real(dp), intent(in) :: stopcr
    real(dp), intent(out) :: sll
    real(dp), intent(inout) :: par(npar), lower(npar), upper(npar)
    interface
       subroutine loglik(model, npar, par, sll)
         import :: dp, argsModel
         implicit none
         integer, intent(in) :: npar
         real(dp), intent(in) :: par(npar)
         type(argsModel), intent(inout) :: model
         real(dp), intent(out) :: sll
       end subroutine loglik
    end interface
    real(dp) :: step(npar)

    conv = 4

    !--------------------------------------------------------------------------------------------
    ! allocating bounds and transforming the parameters to (-infty, infty)
    !--------------------------------------------------------------------------------------------
    call set_bounds(model%bounds, lower, upper, nbd, max(1,npar))
    call transform_par(par, npar, model%bounds,.false.)

    ! step for Nelder-Mead. Based on the transformed parameters
    !step =  minval(0.1d0*abs(par), mask =  abs(par) > 0)
    step = max(0.10d0*abs(par), 0.00025d0)
    where(lower == upper .and. nbd == 2) step = 0.d0

    !-------------------------------------------------------------------------------
    ! Nelder-Mead
    !-------------------------------------------------------------------------------
    call minim(par, step, npar, sll, maxit, iprint, stopcr, loglik, conv, neval, model)

    !-------------------------------------------------------------------------------
    ! transforming the parameters for the original scale
    !------------------------------------------------------------------------------
    call  transform_par(par, npar, model%bounds,.true.)

    return      
  end subroutine optim_nelder

  subroutine optim_lbfgsb(loglik, model, npar, par, nbd, lower, upper,&
       sll, score, iprint, factr, pgtol, maxit, neval, conv, convm)
    !---------------------------------------------------------
    !
    ! Optimization subroutine
    ! The log-likelihood function used here must return the
    ! loglikelihood value and the score vector
    !
    !---------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    integer, intent(in) :: iprint, maxit, npar
    integer, intent(inout) :: nbd(npar)
    integer, intent(out) ::  neval, conv
    real(dp), intent(in) :: factr, pgtol
    real(dp), intent(out) :: sll, score(npar)
    real(dp), intent(inout) :: par(npar), lower(npar), upper(npar)
    character(len = 60), intent(out) :: convm
    interface
       subroutine loglik(model, npar, par, sll, score)
         import :: dp, argsModel
         implicit none
         type(argsModel), intent(inout) :: model
         integer, intent(in) :: npar
         real(dp), intent(in) :: par(npar)
         real(dp), intent(out) :: sll, score(npar)
       end subroutine loglik
    end interface

    integer, parameter ::  mmax = 5 ! default do r
    real(dp) ::  dsave(29)       
    ! wa must have size 2*mn + 11m**2 + 5n + 8m
    real (dp) :: wa(2*mmax*npar+12*mmax*mmax+5*npar+8*mmax)      
    integer :: iwa(3*npar), isave(44)
    character (len=60) :: csave
    logical :: lsave(4)
    integer :: niter

    wa = 0.d0
    iwa = 0
    niter = 0

    ! l-bfgs-b
    !=========================================================
    ! setting lower bounds, upper bounds and fixed values
    !=========================================================
    !     nbd is an integer array of dimension n that must be set by the
    !       user to the type of bounds imposed on the variables:
    !       nbd(i)=0 if x(i) is unbounded,
    !              1 if x(i) has only a lower bound,
    !              2 if x(i) has both lower and upper bounds,
    !              3 if x(i) has only an upper bound.


    convm = 'start'
    model%sco = 0
    call loglik(model, npar, par, sll, score)
    model%sco = 1
    conv = 1

    if(maxit == 0) return

111 continue

    niter = niter + 1

    call setulb(npar,mmax, par,lower,upper,nbd,sll,score,factr,pgtol,wa,iwa,convm,iprint,&
         csave,lsave,isave,dsave)

    if (convm(1:2) .eq. 'fg') then
       call loglik(model, npar, par, sll, score)
       neval = isave(13)
       if(niter > maxit) then
          convm = "max number of iteration reached"
          goto 222
       end if
       goto 111
    endif

    if (convm(1:5) .eq. 'new_x')  goto 111
    !     if task is neither fg nor new_x we terminate execution.
    if(niter > maxit) then
       convm = "max number of iteration reached"
       goto 222
    end if

222 continue

    if(convm(1:4) == 'conv') conv = 0

    return      
  end subroutine optim_lbfgsb


  !*********************************************************************************************
  !
  !                       Subroutines used to allocate vectors and matrices
  !                       related to the models implemented here.    
  !
  !*********************************************************************************************
  subroutine allocate_parvec(vec, length, fix, flags, fval)
    !--------------------------------------------------------------
    !
    ! Helper: subroutine used to allocate polynomials
    !
    ! vec = vector of coefficients
    ! length = size of the vector
    ! fix = number of fixed coefficients
    ! flags = lags of fixed coefficients
    ! fval = fixed values
    !
    !--------------------------------------------------------------
    implicit none
    integer, intent(in) :: fix, length
    integer, intent(in) :: flags(max(1,fix))
    real(dp), intent(in) :: fval(max(1,fix))
    type(vec_parameter), intent(inout) :: vec
    integer :: dummy(max(1,length)), i

    vec%length = length    
    vec%fit =  length - fix
    if(length == 0) return

    call safe_allocate(vec%par, length)
    vec%par = 0.d0
    dummy = 1    

    ! fixed lags and fixed values
    if(fix > 0) then
       dummy(flags) = 0  
       call safe_allocate(vec%flags, fix)
       vec%flags = flags
       vec%par(vec%flags) = fval       
    end if

    ! non-fixed lags 
    if(length - fix > 0) then
       call safe_allocate(vec%lags, length - fix)
       vec%lags = pack((/(i, i = 1,length)/), dummy == 1)
    end if
    return
  end subroutine allocate_parvec

  subroutine allocate_model_ts(model, n, y, gy, inf)
    !---------------------------------------------------------------------------
    !
    ! Helper: subroutine to allocate the time series common to all models
    ! 
    ! each particular model will have extra allocation settings if necessary
    !
    !---------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    integer, intent(in) :: n, inf
    real(dp), intent(in) :: y(n), gy(n)

    model%n = n
    model%inf(1) = inf
    call safe_allocate(model%y, n)
    call safe_allocate(model%gy, n)
    call safe_allocate(model%error,n)      
    model%y = y
    model%gy = gy
    model%error = 0.d0
    return
  end subroutine allocate_model_ts

  subroutine allocate_conditional_ts(cts, n, nreg, xreg, xstart)
    !---------------------------------------------------------------------------
    !
    ! Helper: subroutine to allocate conditional variables
    ! 
    ! each particular model will have extra allocation settings if necessary
    !
    ! ut is the conditional mean or variance
    ! eta = g(ut)
    ! xreg is the corresponding regressor
    !
    !---------------------------------------------------------------------------
    implicit none
    type(conditional_ts), intent(inout) :: cts
    integer, intent(in) :: n, nreg
    real(dp), intent(in) :: xreg(n, max(1,nreg))
    real(dp), intent(in) :: xstart(max(1,nreg))
    cts%nreg = nreg
    call safe_allocate(cts%xstart, max(1,nreg))
    cts%xstart = xstart
    if(nreg > 0) then
       call safe_allocate(cts%xreg, n, nreg)
       cts%xreg = xreg
    end if
    call safe_allocate(cts%ut, n)
    call safe_allocate(cts%eta, n)
    cts%ut = 0.d0    
    cts%eta = 0.d0        
    return
  end subroutine allocate_conditional_ts

  subroutine allocate_model_part(model, fita, alpha,  p, fitar, flagsar, fvar,&
       q, fitma, flagsma, fvma, n, nreg, fitb, flagsb, fvbeta, xreg, xregar, &
       xstart, fitd, d, inf, m, part)
    !---------------------------------------------------------------------------
    !
    ! Helper: subroutine to allocate the base of all models and to set
    !         fixed values 
    ! 
    ! Each particular model will have extra allocation settings if necessary.
    !
    !---------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    integer, intent(in) :: inf, part
    integer, intent(in) :: nreg, p, q, n, m
    integer, intent(in) :: fita, fitb, fitar, fitma, fitd
    integer, intent(in) :: flagsar(max(1,p-fitar)), flagsma(max(1,q-fitma))
    integer, intent(in) :: flagsb(max(1,nreg-fitb)), xregar
    real(dp), intent(in) :: alpha, fvbeta(max(1,nreg-fitb)), xstart(max(1,nreg))
    real(dp), intent(in) :: fvar(max(1,p-fitar)), fvma(max(1,q-fitma)), d
    real(dp), intent(in) :: xreg(n, max(1,nreg))

    !-------------------------------------------------------
    ! prec: if invoked once (part = 1) sets fixprec = 1
    !       in the second call (part = 2) sets fixprec = 0
    !------------------------------------------------------
    model%fixnu = part*(2-part)

    !----------------------------------------------
    ! alpha, beta, ar, ma and d
    !----------------------------------------------
    call allocate_parvec(model%alpha(part), 1, 1-fita, (/1/), (/alpha/))
    call allocate_parvec(model%beta(part), nreg, nreg-fitb, flagsb, fvbeta)
    call allocate_parvec(model%ar(part), p, p-fitar, flagsar, fvar)
    call allocate_parvec(model%ma(part), q, q-fitma, flagsma, fvma)    
    call allocate_parvec(model%d(part), 1, 1-fitd, (/1/), (/d/))

    !---------------------------------------------
    !  If d = 0 uses inf = q
    !---------------------------------------------
    model%inf(part) = max(inf, q)
    if(d == 0.0d0 .and. fitd == 0) model%inf(part) = q

    model%m = m
    !----------------------------------------------
    ! mu, eta and xreg
    !----------------------------------------------
    call allocate_conditional_ts(model%cts(part), n, nreg, xreg, xstart)
    model%cts(part)%xregar = xregar
    return      
  end subroutine allocate_model_part

  subroutine allocate_model1(model, n, y, gy, nreg, xreg, xstart, fitnu, nu, fita, alpha,  &
       fitb, flagsb, fvbeta, xregar, p, fitar, flagsar, fvar, q, fitma, flagsma, fvma, &
       fitd, d, inf, m)    
    !---------------------------------------------------------------------------
    !
    ! Main subroutine to allocate the base of all models with fixed precision
    !
    ! Fixed values are also set at this point
    ! 
    ! each particular model will have extra allocation settings if necessary
    !
    !---------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    integer, intent(in) :: inf
    integer, intent(in) :: fitnu, fita, fitb, fitar, fitma, fitd
    integer, intent(in) :: n, nreg, p, q, m
    integer, intent(in) :: flagsar(max(1,p-fitar)), flagsma(max(1,q-fitma))
    integer, intent(in) :: flagsb(max(1,nreg-fitb)), xregar
    real(dp), intent(in) :: nu, alpha, fvbeta(max(1,nreg-fitb)), d
    real(dp), intent(in) :: fvar(max(1, p-fitar)), fvma(max(1, q - fitma))
    real(dp), intent(in) :: xreg(n, max(1,nreg)), xstart(max(1,nreg))
    real(dp), intent(in) :: y(n), gy(n)
    call allocate_parvec(model%nu, 1, 1-fitnu, (/1/), (/nu/))
    call allocate_model_part(model, fita, alpha,  p, fitar, flagsar, fvar,&
         q, fitma, flagsma, fvma, n, nreg, fitb, flagsb, fvbeta, xreg, xregar,&
         xstart, fitd, d, inf, m, 1)
    call allocate_model_ts(model, n, y, gy, model%inf(1))
    model%npar(1) = fitnu + fita + fitb + fitar + fitma + fitd
    model%npar(2) = 0
    return    
  end subroutine allocate_model1

  subroutine allocate_model2(model, n, y, gy, nreg, xreg1, xreg2, xstart1, xstart2, fita, alpha, &
       fitb, flagb1, flagb2, fvbeta1, fvbeta2, xregar, p, fitar,  flagar1, flagar2, fvar1, fvar2, &
       q, fitma, flagma1, flagma2, fvma1, fvma2, fitd, d, inf, m)

    !---------------------------------------------------------------------------
    !
    ! Main subroutine to allocate the base of all models with fixed precision
    !
    ! Fixed values are also set at this point
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !---------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    integer, intent(in) :: inf(2), m
    integer, intent(in) :: fita(2), fitb(2),  fitar(2), fitma(2), fitd(2)
    integer, intent(in) :: n, nreg(2), p(2), q(2), xregar(2)
    integer, intent(in) :: flagar1(max(1, p(1)-fitar(1)))
    integer, intent(in) :: flagar2(max(1, p(2)-fitar(2)))
    integer, intent(in) :: flagma1(max(1, q(1)-fitma(1)))
    integer, intent(in) :: flagma2(max(1, q(2)-fitma(2))) 
    integer, intent(in) :: flagb1(max(1, nreg(1)-fitb(1)))
    integer, intent(in) :: flagb2(max(1, nreg(2)-fitb(2)))
    real(dp), intent(in) :: alpha(2), d(2)
    real(dp), intent(in) :: fvar1(max(1, p(1)-fitar(1)))
    real(dp), intent(in) :: fvar2(max(1, p(2)-fitar(2)))
    real(dp), intent(in) :: fvma1(max(1, q(1)-fitma(1)))
    real(dp), intent(in) :: fvma2(max(1, q(2)-fitma(2)))    
    real(dp), intent(in) :: fvbeta1(max(1, nreg(1)-fitb(1)))
    real(dp), intent(in) :: fvbeta2(max(1, nreg(2)-fitb(2)))
    real(dp), intent(in) :: xreg1(n, max(1,nreg(1))), xreg2(n, max(1,nreg(2)))
    real(dp), intent(in) :: xstart1(max(1,nreg(1))), xstart2(max(1,nreg(2)))
    real(dp), intent(in) :: y(n), gy(n)
    integer:: infi

    call allocate_parvec(model%nu, 1, 0, (/1/), (/0.d0/))  ! dummy  
    call allocate_model_part(model, fita(1), alpha(1), p(1), fitar(1), flagar1, fvar1, &
         q(1), fitma(1), flagma1, fvma1,  n,  nreg(1), fitb(1), flagb1, fvbeta1, xreg1, xregar(1), &
         xstart1, fitd(1), d(1), inf(1), m, 1)   
    call allocate_model_part(model, fita(2), alpha(2),  p(2), fitar(2), flagar2,  fvar2,&
         q(2), fitma(2), flagma2, fvma2, n, nreg(2), fitb(2), flagb2, fvbeta2, xreg2, xregar(2), &
         xstart2, fitd(2), d(2), inf(2), m, 2)
    infi = max(model%inf(1), model%inf(2))
    call allocate_model_ts(model, n, y, gy,infi)
    model%npar(1) = fita(1) + fitb(1) + fitar(1) + fitma(1) + fitd(1)  
    model%npar(2) = fita(2) + fitb(2) + fitar(2) + fitma(2) + fitd(2)  
    return    
  end subroutine allocate_model2

  subroutine allocate_deta(deta, fita, fitb, fitar, fitma, fitd, n)
    !---------------------------------------------------------------------------
    !
    ! Helper: subroutine to allocate the common parts in the score vector
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !---------------------------------------------------------------------------
    implicit none
    type(deta_d), intent(inout) :: deta
    integer, intent(in) :: n, fita, fitb, fitar, fitma, fitd

    if(fita == 1) call safe_allocate(deta%dalpha, n, 1)
    if(fitb > 0) call safe_allocate(deta%dbeta, n, fitb)
    if(fitar > 0) call safe_allocate(deta%dphi, n, fitar)
    if(fitma > 0) call safe_allocate(deta%dtheta, n, fitma)
    if(fitd == 1) call safe_allocate(deta%dd, n, 1)

    return
  end subroutine allocate_deta

  subroutine allocate_Us(U, fitnu, fita, fitb, fitar, fitma, fitd)
    !---------------------------------------------------------------------------
    !
    ! Helper: subroutine to allocate the common parts in the score vector
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !---------------------------------------------------------------------------
    implicit none
    type(score), intent(inout) :: U
    integer, intent(in) :: fitnu, fita, fitb, fitar, fitma, fitd

    if(fitnu > 0) call safe_allocate(U%Unu,fitnu)
    if(fita == 1) call safe_allocate(U%Ualpha, 1)
    if(fitb > 0) call safe_allocate(U%Ubeta, fitb)
    if(fitar > 0) call safe_allocate(U%Uphi, fitar)
    if(fitma > 0) call safe_allocate(U%Utheta, fitma)
    if(fitd == 1) call safe_allocate(U%Ud, 1)

    return
  end subroutine allocate_Us

  subroutine allocate_SI(model, SI)
    !---------------------------------------------------------------------------
    !
    ! Main subroutine to allocate the common parts in the score vector
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !---------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(in) :: model
    type(argsSI), intent(inout) :: SI

    ! T1 and h1
    call safe_allocate(SI%T(1)%par, model%n)
    call safe_allocate(SI%h(1)%par, model%n)

    ! Urho
    if(model%sco == 1) call allocate_Us(SI%U(1), model%nu%fit, model%alpha(1)%fit,&
         model%beta(1)%fit, model%ar(1)%fit, model%ma(1)%fit, model%d(1)%fit)

    ! deta1_drho
    call allocate_deta(SI%deta(1,1), model%alpha(1)%fit,&
         model%beta(1)%fit, model%ar(1)%fit, model%ma(1)%fit, model%d(1)%fit, model%n)

    if(model%fixnu == 1) return    

    ! T2 and h2
    call safe_allocate(SI%T(2)%par, model%n)
    call safe_allocate(SI%h(2)%par, model%n)

    ! Ulambda
    if(model%sco == 1) call allocate_Us(SI%U(2), 0, model%alpha(2)%fit,&
         model%beta(2)%fit, model%ar(2)%fit, model%ma(2)%fit, model%d(2)%fit)

    ! deta1_dlambda
    call allocate_deta(SI%deta(1,2), model%alpha(2)%fit,&
         model%beta(2)%fit, model%ar(2)%fit, model%ma(2)%fit, model%d(2)%fit, model%n)
    ! deta2_drho
    call allocate_deta(SI%deta(2,1), model%alpha(1)%fit,&
         model%beta(1)%fit, model%ar(1)%fit, model%ma(1)%fit, model%d(1)%fit, model%n)
    ! deta2_dlambda
    call allocate_deta(SI%deta(2,2),model%alpha(2)%fit,&
         model%beta(2)%fit, model%ar(2)%fit, model%ma(2)%fit, model%d(2)%fit, model%n)

    return
  end subroutine allocate_SI

  !******************************************************************
  !
  !  Subroutines used to initialize/update the vector of parameters  
  !
  !******************************************************************
  subroutine start_par1(par, model, part)
    !-----------------------------------------------------------------------------------
    !
    ! Helper: Initialize the non-fixed parameters
    !
    ! This subroutine is called to set the parameter values during the log-likelihood
    ! evaluation (during and after estimation process)
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !------------------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    real(dp), intent(in) :: par(sum(model%npar))
    integer, intent(in) :: part
    integer :: n1, n2

    n1 = (2 - part) + (part - 1)*(model%npar(1) + 1)
    n2 = n1 - 1

    ! alpha
    n1 = n2 + 1
    n2 = n2 + model%alpha(part)%fit
    if(n2 >= n1) model%alpha(part)%par = par(n1:n2)

    ! beta
    n1 = n2 + 1
    n2 = n2 + model%beta(part)%fit
    if(n2 >= n1) model%beta(part)%par(model%beta(part)%lags) = par(n1:n2)

    ! phi-ar
    n1 = n2 + 1
    n2 = n2 + model%ar(part)%fit
    if(n2 >= n1) model%ar(part)%par(model%ar(part)%lags) = par(n1:n2)

    ! theta-ma
    n1 = n2 + 1
    n2 = n2 + model%ma(part)%fit
    if(n2 >= n1) model%ma(part)%par(model%ma(part)%lags) = par(n1:n2)

    ! d
    n1 = n2 + 1
    n2 = n2 + model%d(part)%fit
    if(n2 >= n1) model%d(part)%par = par(n1:n2)

    ! for compatibility with models where the dispersion depends on t,
    ! nu is always the last parameter
    if(model%fixnu == 1 .and. part == 1) then
       ! nu
       if(model%nu%fit > 0) then
          n1 = model%npar(1) - model%nu%fit + 1
          n2 = model%npar(1)
          model%nu%par = par(n1:n2)
       end if
    end if

    return   
  end subroutine start_par1

  subroutine start_par2(par, model, vc, part)
    !-----------------------------------------------------------------------------------
    !
    ! Helper: Initialize the non-fixed parameters
    !
    ! This subroutine is called to set the parameter values during the log-likelihood
    ! evaluation (during and after estimation process)
    !
    ! each particular model will have extra allocation settings if necessary
    !
    !------------------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    real(dp), intent(in) :: par(sum(model%npar))
    integer, intent(in) :: part
    real(dp), intent(out) :: vc(0:model%inf(part))
    real(dp) :: theta(0:model%ma(part)%length)

    call start_par1(par, model, part)

    theta(0) = 1.d0
    if(model%ma(part)%length > 0) theta(1:model%ma(part)%length) = model%ma(part)%par
    call vc_f(model%d(part)%par(1), theta, model%ma(part)%length, model%inf(part), vc)

    return
  end subroutine start_par2

  !**************************************************************************************************
  !
  !                  Subroutines used to process information in and out the model
  !
  !**************************************************************************************************
  subroutine get_model(model, n, y, ylower, yupper, ystart, gy, nreg, xreg, xstart, &
       link, escale, skippar, npar, par, fixa, alpha, fixB, flagsb, fvbeta, &
       p, fixar, flagsar, fvar, xregar, q, fixma, flagsma, fvma, fixd, d, fixnu, nu, &
       inf, llk, sco, info, m)
    !-----------------------------------------------------------------------------------------
    !
    !  Subrotuine used to pass the values enter in the main program to the user defined
    !  variable that will be passed to the  generic subroutines
    !
    !-----------------------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    integer, intent(in) :: n, nreg, link(2), inf, llk, sco, info, m, skippar, escale
    integer, intent(in) ::  npar,  fixa, fixB, fixar, fixma, fixd, fixnu
    integer, intent(in) :: p, q, xregar
    integer, intent(in) :: flagsb(max(1,fixB)),flagsar(max(1,fixar)), flagsma(max(1,fixma))
    real(dp), intent(in) :: y(n), xreg(n, max(1, nreg)), ylower, yupper, ystart
    real(dp), intent(in) :: par(npar), alpha, d, nu, fvbeta(max(1,fixB))
    real(dp), intent(in):: fvar(max(1,fixar)), fvma(max(1,fixma)), xstart(max(1, nreg))
    real(dp), intent(inout) :: gy(n)
    integer :: t

    model%m = m
    model%ystart = ystart
    model%escale = escale

    !-------------------------------------------------------------
    ! log-likelihood, score vector and information matrix
    !-------------------------------------------------------------
    model%llk = llk
    model%sco = sco
    model%info = info

    !----------------------------------------------------------------
    ! setting the link parameters and calculating g(y)
    !----------------------------------------------------------------
    if(.not. allocated(model%argsL)) allocate(model%argsL(2))        
    model%argsL(1:2)%link = link
    model%argsL(1:2)%lower = ylower
    model%argsL(1:2)%upper = yupper         
    do t = 1, n
       gy(t) = linkfun(y(t), model%argsL(2))
    end do

    !-------------------------------------------------------------
    ! allocating the time series and parameters
    ! setting fixed values/lags of parameters
    !-------------------------------------------------------------
    call allocate_model(model, n, y, gy, nreg, xreg, xstart,1-fixnu, nu, 1-fixa, alpha, &
         nreg-fixb, flagsb, fvbeta, xregar, p, p-fixar, flagsar, fvar, q, q-fixma, flagsma, fvma,&
         1-fixd, d, inf, m)

    !--------------------------------------------------------------
    !  Setting the parameter's non-fixed values:
    !  alpha, beta, ar, ma and d
    !  BARC models must skip this step an call start_par later
    !--------------------------------------------------------------
    if(skippar == 0) call start_par(par, model, 1)

    if(sco + info == 0 ) return
    !---------------------------------------------------------------------
    ! allocating score-vector, Information matrix and related matrices
    !---------------------------------------------------------------------
    call allocate_SI(model, model%SI)

    return    
  end subroutine get_model

  subroutine return_model(model, n, mu, eta, error, inf, extra, nd, D, T, ne, E, h)
    implicit none
    type(argsModel), intent(in) :: model
    integer, intent(in) :: n, nd, ne, extra
    real(dp), intent(out), dimension(n) :: mu, eta, error
    integer, intent(out) :: inf
    real(dp), intent(out) :: T(max(n*extra, 1)), h(max(n*extra, 1))
    real(dp), intent(out) :: D(max(n*extra, 1), max(nd*extra,1))
    real(dp), intent(out) :: E(max(n*extra, 1), max(ne*extra, 1))

    mu = model%cts(1)%ut
    eta = model%cts(1)%eta
    error = model%error
    inf = model%inf(1)

    D = 0.d0
    E = 0.d0
    h = 0.d0
    if(extra == 1 .and. (model%sco + model%info) > 0) then       
       call fill_D(model%SI, model%alpha(1)%fit, model%beta(1)%fit, model%ar(1)%fit, &
            model%ma(1)%fit, model%d(1)%fit, n, nd, D)
       T = model%SI%T(1)%par
       if(model%info == 1 .and. (model%npar(1) - model%nu%fit) > 0) E = model%SI%E
       if(model%sco == 1 .and. (model%npar(1) - model%nu%fit) > 0) h = model%SI%h(1)%par
    end if
    return
  end subroutine return_model

  subroutine final_model(model, npar, par, fixnu, n, mu, eta, error, nnew, nreg, &   
       xnew, ynew, inf, sll, sco, U, info, K, extra, Drho, T, E, h, &        
       llk_dist, dllk_dist, Ed2llk_dist)
    !-----------------------------------------------------------------------------
    !
    !   Reports the final model and related matrices
    !
    !------------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    integer, intent(in) :: npar, fixnu, n, nnew, nreg, info, extra, sco
    integer, intent(out) :: inf
    real(dp), intent(in) :: par(npar), xnew(max(1,nnew), max(1, nreg))
    real(dp), intent(out) :: mu(n), eta(n), error(n), ynew(max(1,nnew))
    real(dp), intent(out) :: sll, U(max(1, npar*sco))
    real(dp), intent(out) :: K(max(1, npar*info), max(1, npar*info))
    real(dp), intent(out) :: T(max(1, n*extra)), h(max(1, n*extra))
    real(dp), intent(out) :: Drho(max(1,n*extra), max(1,npar - 1 + fixnu))
    real(dp), intent(out) :: E(max(1,n*extra), 1+2*(1-fixnu)*extra)
    !-----------------------------------------------------------------------
    ! llk_dist and dllk_dist: arguments in loglik_generic
    ! Ed2llk_dist: arguments in K_generic
    !-----------------------------------------------------------------------         
    interface
       function llk_dist(m, n, y, mu, nu, argsD) result(sll)
         import :: dp, argsDist
         implicit none
         integer, intent(in) :: m, n
         real(dp), intent(in) :: y(n), mu(n), nu
         type(argsDist), intent(in) :: argsD
         real(dp) :: sll
       end function llk_dist
    end interface
    interface
       subroutine dllk_dist(m,n,y,n1,mut,skipmu,n2,nut,skipnu,dllmu,dllnu,argsD)
         import :: dp, argsDist
         implicit none
         type(argsDist), intent(in) :: argsD
         integer, intent(in) :: m, n, n1, n2, skipmu, skipnu
         real(dp), intent(in) :: y(n), mut(n1), nut(n2)
         real(dp), intent(out) :: dllmu(n*(1-skipmu)+skipmu)
         real(dp), intent(out) :: dllnu(n*(1-skipnu)+skipnu)
       end subroutine dllk_dist
    end interface
    interface
       subroutine Ed2llk_dist(m,n,n1,mut,i1,n2,nut,i2,E,argsD)
         import :: dp, argsDist
         implicit none
         integer, intent(in) :: m, n, n1, n2, i1, i2
         real(dp), intent(in) :: mut(n1), nut(n2)
         real(dp), intent(out) :: E(n,max(1,i1+i2+i1*i2))
         type(argsDist), intent(in) :: argsD
       end subroutine Ed2llk_dist
    end interface
    real(dp), allocatable :: vc(:)
    real(dp) :: Uaux(npar)

    !------------------------------------------------------------------------------
    ! reporting score vector and information matrix (if required)
    !------------------------------------------------------------------------------
    model%llk = 1
    if(sco + info > 0) then
       ! setting model%sco  = 1 so D and T will be calculated
       model%sco = max(sco,info)
       model%info = info
       call allocate_SI(model, model%SI)
    end if

    sll = 0.d0
    U = 0.d0
    K = 0.d0
    Uaux = 0.d0      
    !---------------------------------------------------------------------------------
    ! Reports the final model.
    ! Calculates: mu, eta, rt, -sll (if llk = 1) and -U (if sco = 1).
    !---------------------------------------------------------------------------------
    call loglik_generic(llk_dist, dllk_dist, model, npar, par, sll, Uaux)      

    if(sco == 1) U = Uaux      

    if(info == 0)  goto 100
    !------------------------------------------------------------------------------
    ! Information matrix
    !------------------------------------------------------------------------------
    call K_generic(Ed2llk_dist, model%SI, model%cts(1)%ut, model%alpha(1)%fit, &
         model%beta(1)%fit, model%ar(1)%fit, model%ma(1)%fit, model%d(1)%fit, model%nu%fit,&
         model%nu%par(1), model%m, model%n, npar, K, model%argsD)

100 continue
    !---------------------------------------------------------------------------------
    ! calling the subroutine to fill the matrices with the calculated values
    !--------------------------------------------------------------------------------
    call return_model(model, n, mu, eta, error, inf, extra, size(Drho,2), Drho,&
         T, size(E,2), E, h)

    !---------------------------------------------------------------------------------  
    ! positive log-likelihood and score vector
    !--------------------------------------------------------------------------------- 
    sll = -sll
    U = -U

    if(nnew == 0) return
    !---------------------------------------------------------------------------------
    !  out-of-sample forecast 
    !---------------------------------------------------------------------------------    
    call safe_allocate(vc, 0, model%inf(1))
    call start_par(par, model, vc, 1)
    call mu_forecast(model, vc, nnew, xnew, ynew)

    return
  end subroutine final_model

  !**************************************************************************************************
  !
  !                  Subroutines used to simulate and predict 
  !
  !**************************************************************************************************
  subroutine sim_model(rdist, n, burn, np, pdist, alpha, nreg, beta, p, phi, &
       q, theta, d, linkg, xreg, xregar, yt, ylower, yupper, ystart, xstart,&
       mu, eta, error, escale, ns, seed, rngtype, inf, rev, argsL)
    !#---------------------------------------------------------------------
    !#
    !#  Simulating a F-ARFIMA model  (F = Beta, Gamma, Kumaraswamy)
    !#
    !#---------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n, burn, nreg, p, q, linkg(2), ns, np
    integer, intent(inout) :: seed(ns), rngtype, xregar, escale, inf
    integer, intent(inout) :: rev
    real(dp), intent(in) :: alpha, pdist(np), d, ylower, yupper
    real(dp), intent(inout) :: beta(max(1,nreg)), phi(max(1,p)), theta(max(1,q))
    real(dp), intent(inout) :: xreg(n+burn,max(1, nreg))
    real(dp), intent(inout) :: ystart, xstart(max(1,nreg))
    real(dp), intent(inout) :: yt(n+burn), mu(n+burn), eta(n+burn), error(n+burn)
    type(argslink), intent(inout) :: argsL(2)
    interface
       function rdist(npar, par, rng) result(fnval)
         import :: dp, rng_t
         implicit none
         integer, intent(in) :: npar
         real(dp), intent(in) :: par(npar)
         type(rng_t), intent(inout) :: rng
         real(dp) ::  fnval
       end function rdist
    end interface
    integer ::  t, i
    type(rng_t) :: rng
    real(dp) :: vc(0:inf), ytemp, xb, gytemp, gy(n+burn)

    !---------------------------------------------
    !  If d = 0 uses inf = q
    !---------------------------------------------
    inf = max(inf, q)
    if(d == 0.0d0) inf = q

    ! revision required
    rev = 1

    !# check for regressors
    if(nreg == 0) then
       xreg = 0.d0
       beta = 0.d0
    end if

    !# check for AR component
    if(p == 0) phi = 0.d0

    !# check for MA component
    if(q == 0) theta = 0.d0

    !# initializing variables in the model
    eta = 0.d0
    mu = 0.d0
    yt = 0.d0
    gy = 0.d0
    error = 0.d0

    argsL(1:2)%link = linkg ! g1(mu) and g2(y)
    argsL(1:2)%lower = ylower
    argsL(1:2)%upper = yupper   

    call rng_seed(rng, ns, seed, rngtype)

    call vc_f(d, (/1.d0,theta/), q, inf, vc)

    ! starting values 
    ytemp = 0.d0
    xb = 0.d0
    if(p > 1) then
       if(ylower < ystart .and. ystart < yupper) ytemp = linkfun(ystart, argsL(2))
       if(xregar == 1) xb = sum(xstart*beta)
    end if

    do t = 1, (n+burn)      
       !#------------------------------------
       !# eta(t) = alpha + x*b + AR + MA
       !#------------------------------------
       eta(t) = alpha + sum(xreg(t,:)*beta)

       do i = 1,p          
          ! updating the auxiliar variables ytemp and xb
          if(t-i > 0) then
             ! update g2(y)
             ytemp = gy(t-i)
             ! check if xreg is used in AR recursion
             if(xregar == 1) xb = sum(xreg(t-i,:)*beta)
          end if
          eta(t) = eta(t) + (ytemp - xb)*phi(i)
       end do

       do i = 1, min(t-1, inf)
          ! sum(c(j)*r(t-j)), only exists if d > 0 or q > 0          
          ! if inf < 1 does not enter the loop
          eta(t) =  eta(t) + vc(i)*error(t-i) 
       end do

       !#---------------------------
       !# mu(t) = g^{-1}(eta(t))
       !#---------------------------
       mu(t) = linkinv(eta(t), argsL(1))
       if(mu(t) < ylower .or. mu(t) > yupper) return

       ! to avoid the bounds:
       if(mu(t) <= ylower) then
          mu(t) = ylower + epsilon(1.d0)
          eta(t) = linkfun(mu(t), argsL(1))
       elseif(mu(t) >= yupper) then
          mu(t) = yupper - epsilon(1.d0)
          eta(t) = linkfun(mu(t), argsL(1))
       end if

       !#------------------------
       !# y(t) ~ f(y, mu, pdist)
       !#------------------------
       yt(t) = rdist(np+1, (/mu(t), pdist/), rng)
       if(yt(t) <= ylower) then
          yt(t) = ylower + epsilon(1.d0)
       elseif(yt(t) >= yupper) then
          yt(t) = yupper - epsilon(1.d0)
       end if
       gy(t) =  linkfun(yt(t), argsL(2))

       if(escale == 0) then
          error(t) = yt(t) - mu(t)
       else
          if(argsL(1)%link /= argsL(2)%link) then
             gytemp = linkfun(yt(t), argsL(1))
          else
             gytemp = gy(t)
          end if
          error(t) = gytemp - eta(t)
       end if
    end do

    rev = 0
    return    
  end subroutine sim_model


  !*********************************************************************************************
  !
  !                    Distribution related subroutines
  !
  !*********************************************************************************************
  function llk_beta(m, n, y, mu, nu, argsD) result(sll)
    !-----------------------------------------------------------------------
    !
    ! Computes 
    !           sll = sum(log(f(y| mu, nu)))
    !  for beta regression models, where f(y| mu, nu) is the beta density
    !
    !-------------------------------------------------------------------------
    ! Input:
    !   m = starting time to calculate the log-likelihood
    !       for t < m + 1, llk = 0
    !
    !   n = sample size
    !
    !   y = observed time series
    !
    !   mu = mean (same size of y)
    !
    !   nu = dispersion parameter (assumed to be constant)
    !
    !   argsD = dummy argument used for compatibility with other
    !           distributions
    !
    ! Output:
    !   sll = sum of the log-likelihood from times m to n. 
    !-----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n
    real(dp), intent(in) :: y(n), mu(n), nu
    type(argsDist), intent(in) :: argsD  ! dummy, for compatibility
    real(dp) :: sll
    integer :: t

    sll = argsD%arg1! Dummy to avoid compiler warning
    sll = 0.d0
    do t = (m +1),n
       sll = sll + dbeta(y(t), 2, (/mu(t), nu/), .true.)
    end do
    return
  end function llk_beta

  subroutine dllk_beta(m,n,y,n1,mut,skipmu,n2,nut,skipnu,dllmu,dllnu,argsD)
    !-------------------------------------------------------------------------
    !
    !  Calculates dl/dmu and dl/dnu for beta regression models
    !
    !  This subroutine works for constant and for time varying parameters.
    !
    !-------------------------------------------------------------------------
    ! Input:
    !   m = starting time to calculate the derivatives
    !       for t < m, dl/dmu = dl/dnu = 0
    !
    !   n = sample size
    !
    !   y = observed time series
    !
    !   n1 = 1, if mu is fixed and n1 = n if mu is time varying
    !
    !   mut = mean
    !
    !   skipmu = indicates if dl/dmu must be calculated
    !            skipmu = 1 indicates that mu is fixed in the model
    !
    !   n2 = 1, if mu is fixed and n2 = n if nu is time varying
    !
    !   nut = dispersion parameter
    !
    !   skipnu = indicates if dl/dnu must be calculated
    !            skipnu = 1 indicates that nu is fixed in the model
    !
    !   argsD = dummy argument used for compatibility with other
    !           distributions
    !
    ! Output:
    !   dllmu and dldnu = the derivatives of the log-likelihood function
    !                     with respect to mu and nu, respectively
    !-------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n, n1, n2, skipmu, skipnu
    real(dp), intent(in) :: y(n), mut(n1), nut(n2)
    type(argsDist), intent(in) :: argsD  ! dummy, for compatibility      
    real(dp), intent(out) :: dllmu(n*(1-skipmu)+skipmu)
    real(dp), intent(out) :: dllnu(n*(1-skipnu)+skipnu)
    real(dp) :: ystar(n), mustar(n), nu(n), mu(n), dgn(n)
    integer :: t

    dllmu = argsD%arg1! Dummy to avoid compiler warning
    dllmu = 0.d0
    dllnu = 0.d0
    if(n1 + n2 == 0) return

    mu = mut(n1)
    nu = nut(n2)
    if(n1 > 1) mu = mut
    if(n2 > 1) nu = nut

    !-------------------------------------------------------
    ! ystar and mustar
    !-------------------------------------------------------
    ystar = log(y) - log(1.d0-y)
    mustar = digamma(mu(1)*nu(1)) - digamma((1.d0 - mu(1))*nu(1))
    if(n1 > 1 .or. n2 > 1) then
       do t = (m+1), n
          mustar(t) = digamma(mu(t)*nu(t)) - digamma((1.d0 - mu(t))*nu(t))
       end do
    end if

    !-------------------------------------------------------
    ! dl/dnu
    !-------------------------------------------------------
    if(skipnu == 1) goto 200
    dgn = digamma(nu(1))
    do t = (m +1), n
       if(n2 > 1) dgn(t) = digamma(nu(t))
       dllnu(t) = mu(t)*(ystar(t)-mustar(t)) + log(1.d0-y(t)) - &
            digamma((1.d0 - mu(t))*nu(t)) + dgn(t)
    end do

200 continue
    !-------------------------------------------------------------
    ! dl/dmu
    !--------------------------------------------------------------
    if(skipmu == 1) goto 300
    do t = (m +1), n
       dllmu(t) = nu(t)*(ystar(t) - mustar(t))
    end do

300 continue
    return
  end subroutine dllk_beta

  subroutine Ed2llk_beta(m,n,n1,mut,i1,n2,nut,i2,E,argsD)
    !-------------------------------------------------------------------------
    !
    !  Calculates
    !           E(dl2/dmu2), E(dl2/dmudnu) and E(dl2/dnu2) 
    !  the exepcted values of the second derivatives of the log-likelihood 
    !  for beta regression models
    !
    !-------------------------------------------------------------------------
    ! Input:
    !   m = starting time to calculate the loglikelihood
    !       for t < m, sll = 0
    !
    !   n = sample size
    !
    !   n1 = size of mut (n1 = 1 when mu is constant)
    !
    !   mu = mean
    !
    !   i1 = 1, if mu is time varying or a constant to be estimated
    !
    !   n2 = size of nut (n2 = 1 when nu is constant)
    !
    !   nu = dispersion parameter
    !
    !   i2 = 1, if nu is time varying or a constant to be estimated
    !
    !   argsD = dummy argument used for compatibility with other
    !           distributions
    !
    ! Output:
    !    E = matrix with the expected values
    !-------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n, n1, n2, i1, i2
    real(dp), intent(in) :: mut(n1), nut(n2)
    real(dp), intent(out) :: E(n,max(1,i1+i2+i1*i2))
    type(argsDist), intent(in) :: argsD  ! dummy, for compatibility      
    real(dp) ::  nu(n), mu(n), psi, psi1, psi2
    integer :: t, ifault

    E = 0.d0
    if(i1 + i2 == 0) return

    if(n1*n2 == 1) then
       !---------------------------------------
       !    mu and nu don't vary with time
       !---------------------------------------
       psi1 = trigamma(mut(1)*nut(1),ifault)
       psi2 = trigamma((1.d0 - mut(1))*nut(1), ifault)       
       if(i1 == 1) E(:,i1) = (psi1 + psi2)*nut(1)**2
       if(i1*i2 == 1) E(:,i1+i1*i2) = ((psi1+psi2)*mut(1) - psi2)*nut(1)
       if(i2 == 1) then
          psi = trigamma(nut(1), ifault)
          E(:,i1+i2+i1*i2) =  psi1*mut(1)**2 + psi2*(1.d0 - mut(1))**2 - psi
       end if
       return
    end if
    !---------------------------------------
    !    mu and/or nu vary with time
    !---------------------------------------
    psi = 0.d0
    mu = mut(n1)
    nu = nut(n2)
    if(n1 > 1) mu = mut
    if(n2 > 1) nu = nut
    if(n2*i2 == 1) psi =  trigamma(nu(1), ifault) 
    do t = (m+1),n     
       psi1 = trigamma(mu(t)*nu(t),ifault)
       psi2 = trigamma((1.d0 - mu(t))*nu(t), ifault)       
       if(i1 == 1) E(t,i1) = (psi1 + psi2)*nu(t)**2
       if(i1*i2 == 1) E(t,i1+i1*i2) = ((psi1+psi2)*mu(t) - psi2)*nu(t)
       if(i2 == 1) then
          if(n2 > 1) psi =  trigamma(nu(t), ifault)           
          E(t,i1+i2+i1*i2) =  psi1*mu(t)**2 + psi2*(1.d0 - mu(t))**2 - psi
       end if
    end do

    return    
  end subroutine Ed2llk_beta
  !--------------------------- BETA (END) --------------------


  function llk_kuma(m, n, y, mu, nu, argsD) result(sll)
    !-----------------------------------------------------------------------
    ! Computes the sum(log(f(y| mu, nu))) for Kumaraswamy regression models
    !-----------------------------------------------------------------------
    implicit none    
    integer, intent(in) :: m, n
    real(dp), intent(in) :: y(n), mu(n), nu
    type(argsDist), intent(in) :: argsD
    real(dp) :: sll
    integer :: t
    !-----------------------------------------------------
    ! log-likelihood for Kumaraswamy regression models
    !-----------------------------------------------------
    sll = 0.d0
    do t = (m +1),n
       sll = sll + dkuma(y(t), 5, (/mu(t), nu, argsD%arg1, argsD%lower, argsD%upper/), .true.)
    end do
    return
  end function llk_kuma

  subroutine dllk_kuma(m,n,y,n1,mut,skipmu,n2,nut,skipnu,dllmu,dllnu, argsD)
    !-------------------------------------------------------------------------
    !
    !  Calculates dl/dmu and dl/dnu for Kumaraswamy regression models
    !
    !   n1 = 1, if mu is fixed and n1 = n if mu is time varying
    !   n2 = 1, if nu is fixed and n2 = n if nu is time varying    
    !
    !-------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n, n1, n2, skipmu, skipnu
    real(dp), intent(in) :: y(n), mut(n1), nut(n2)
    real(dp), intent(out) :: dllmu(n*(1-skipmu)+skipmu)
    real(dp), intent(out) :: dllnu(n*(1-skipnu)+skipnu)
    type(argsDist), intent(in) :: argsD
    real(dp) :: ct(n), mu1, y01(n)
    real(dp) :: delta(n), nu(n), mu01(n)
    integer :: t

    dllmu = 0.d0
    dllnu = 0.d0
    if(n1 + n2 == 0) return

    y01 = (y - argsD%lower) / (argsD%upper - argsD%lower)    
    mu01 = (mut(n1) - argsD%lower) / (argsD%upper - argsD%lower)
    nu = nut(n2)
    if(n1 > 1) mu01 = (mut - argsD%lower) / (argsD%upper - argsD%lower)
    if(n2 > 1) nu = nut      

    delta = 0.d0
    ct = 0.d0
    !-----------------------------------------
    ! calculating delta(t) and c(t)
    !-----------------------------------------
    if(n1*n2 == 1) then
       !---------------------------------------------
       ! mu and nu are both fixed values
       !---------------------------------------------      
       mu1 =  1.d0 - mu01(1)**nu(1)
       delta = log(1.d0 - argsD%arg1) / log(mu1)
       ct = mu01(1)**(nu(1)-1.d0) / (mu1*log(mu1))
       ct = ct * (delta * log(1.d0 - y01**nu(1)) + 1.d0)
    else     
       !---------------------------------------------
       ! at least one varies with time
       !---------------------------------------------
       do t = (m + 1), n
          mu1 =  1.d0 - mu01(t)**nu(t)
          delta(t) = log(1.d0 - argsD%arg1) / log(mu1)       
          ct(t) = mu01(t)**(nu(t) - 1.d0) / (mu1 * log(mu1))
          ct(t) = ct(t) * (delta(t) * log(1.d0 - y01(t)**nu(t)) + 1.d0)
       end do
    end if

    !-------------------------------------------------------
    ! dl/dnu
    !-------------------------------------------------------
    if(skipnu == 1) goto 200
    do t = (m + 1), n
       dllnu(t) =  1/nu(t) + log(y01(t)) + ct(t)*mu01(t)*log(mu01(t)) - &
            (delta(t) - 1.0d0) * y01(t)**nu(t) * log(y01(t))/(1.d0 - y01(t)**nu(t)) 
    end do

200 continue
    !-------------------------------------------------------------
    ! dl/dmu
    !--------------------------------------------------------------
    if(skipmu == 1) goto 300
    dllmu = nu * ct /  (argsD%upper - argsD%lower)

300 continue
    return
  end subroutine dllk_kuma

  subroutine Ed2llk_kuma(m,n,n1,mut,i1,n2,nut,i2,E, argsD)
    !-------------------------------------------------------------------------
    !
    !  Calculates
    !           E(dl2/dmu2), E(dl2/dmudnu) and E(dl2/dnu2) 
    !  for Kumaraswamy regression models
    !
    !   n1 = 1, if mu is fixed and n1 = n if mu is time varying
    !   n2 = 1, if mu is fixed and n2 = n if nu is time varying   
    !   i1/i2 = 0, if mu/nu are fixed values provided by the user
    !
    !-------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n, n1, n2, i1, i2
    real(dp), intent(in) :: mut(n1), nut(n2)
    real(dp), intent(out) :: E(n,max(1,i1+i2+i1*i2))
    type(argsDist), intent(in) :: argsD    
    real(dp) :: mu1, tau1(n), tau2(n), aux
    real(dp) :: delta(n), nu(n), mu01(n), dgd, ddgd, s
    integer :: t, ifault
    real(dp), parameter :: em = 0.57721566490153286061d0
    real(dp), parameter :: k0 = 0.82368066085287928d0

    E = 0.d0
    if(i1 + i2 == 0) return

    mu01 = (mut(n1) - argsD%lower) / (argsD%upper - argsD%lower)
    nu = nut(n2)
    if(n1 > 1) mu01 = (mut - argsD%lower) / (argsD%upper - argsD%lower)
    if(n2 > 1) nu = nut

    delta = 0.d0
    tau1 = 0.d0
    tau2 = 0.d0
    aux = 0.d0
    if(n1*n2 == 1) then
       !---------------------------------------------
       ! mu and nu are both fixed values
       !---------------------------------------------      
       mu1 =  1.d0 - mu01(1)**nu(1)
       delta = log(1.d0 - argsD%arg1) / log(mu1)       
       tau1 = mu01(1)**(nu(1) - 2.d0)/(mu1*log(mu1))
       tau2 = tau1**2*mu01(1)**2
    else
       !---------------------------------------------
       ! at least one varies with time
       !---------------------------------------------             
       do t = (m + 1), n
          mu1 =  1.d0 - mu01(t)**nu(t)
          delta(t) = log(1.d0 - argsD%arg1) / log(mu1)       
          tau1(t) = mu01(t)**(nu(t) - 2.d0)/(mu1*log(mu1))
          tau2(t) = tau1(t)**2*mu01(t)**2
       end do
    end if

    if(n1*n2 == 1) then
       !---------------------------------------
       !    mu and nu don't vary with time
       !---------------------------------------
       if(i1 == 1) E(:,i1) = tau2(m+1)*(nu(1)/(argsD%upper - argsD%lower))**2
       ! setting the values used by dl2/dmudnu and dl2/dnu2
       if(i2 == 1) then
          dgd = digamma(delta(m+1))
          ddgd = trigamma(delta(m+1), ifault)          
          mu1 = mu01(1)*log(mu01(1))
          ! digamma(x + 1) = digamma(x) + 1/x
          aux = (1.d0 - dgd - 1.d0/delta(m+1) - em) / (delta(m+1) - 1.d0)
          aux = delta(m+1)*tau1(m+1)*aux
          if(i1 == 1) then
             E(:,i1+i1*i2) = E(:,i1)*mu1*(argsD%upper - argsD%lower)/nu(1)
             E(:,i1+i1*i2) = E(:,i1+i1*i2) + aux*mu01(1)/(argsD%upper - argsD%lower)
          end if
          s = 2 * aux * mu1 * mu01(1) / nu(1)      
          E(:,i1+i2+i1*i2) = 1.d0/nu(1)**2 + tau2(m+1)*mu1**2  + s         
          s = (dgd*(dgd + 2.d0*(em - 1.d0)) - ddgd + k0)/nu(1)**2
          s = s*delta(m+1)/(delta(m+1) - 2.d0)
          E(:,i1+i2+i1*i2) = E(:,i1+i2+i1*i2) + s
       end if
       return
    end if

    !---------------------------------------
    !    mu and/or nu vary with time
    !---------------------------------------
    do t = (m+1),n        
       if(i1 == 1) E(t,i1) = tau2(t)*(nu(t)/(argsD%upper - argsD%lower))**2
       ! setting the values used by dl2/dmudnu and dl2/dnu2
       if(i2 == 1) then          
          dgd = digamma(delta(t))
          ddgd = trigamma(delta(t), ifault)               
          mu1 = mu01(t)*log(mu01(t))
          ! digamma(x + 1) = digamma(x) + 1/x
          aux = (1.d0 - dgd - 1.d0/delta(t) - em) / (delta(t) - 1.d0)
          aux = delta(t)*tau1(t)*aux
          if(i1 == 1) then
             E(t,i1+i1*i2) = E(t,i1)*mu1*(argsD%upper - argsD%lower)/nu(t)
             E(t,i1+i1*i2) = E(t,i1+i1*i2) + aux*mu01(t)/(argsD%upper - argsD%lower)
          end if
          s = 2 * aux * mu1 * mu01(t) / nu(t)      
          E(t,i1+i2+i1*i2) = 1.d0/nu(t)**2 + tau2(t)*mu1**2  + s         
          s = (dgd*(dgd + 2.d0*(em - 1.d0)) - ddgd + k0)/nu(t)**2
          s = s*delta(t)/(delta(t) - 2.d0)
          E(t,i1+i2+i1*i2) = E(t,i1+i2+i1*i2) + s                           
       end if
    end do

    return    
  end subroutine Ed2llk_kuma
  !--------------------------- Kumaraswamy (END) --------------------

  function llk_gamma(m, n, y, mu, nu, argsD) result(sll)
    !-----------------------------------------------------------------------
    !
    ! Computes 
    !           sll = sum(log(f(y| mu, nu)))
    ! for gamma regression models, where f(y| mu, nu) is the gamma density
    !
    !-------------------------------------------------------------------------
    ! Input:
    !   m = starting time to calculate the loglikelihood
    !       for t < m, sll = 0
    !
    !   n = sample size
    !
    !   y = observed time series
    !
    !   mu = mean (same size of y)
    !
    !   nu = dispersion parameter (assumed to be constant)
    !
    ! Output:
    !   sll = sum of the log-likelihood from times m to n. 
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Computes the sum(log(f(y| mu, nu))) for gamma regression models
    !-----------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n
    real(dp), intent(in) :: y(n), mu(n), nu
    type(argsDist), intent(in) :: argsD  ! dummy, for compatibility      
    real(dp) :: sll
    integer :: t

    sll = argsD%arg1! Dummy to avoid compiler warning
    !-----------------------------------------
    ! log-likelihood for gamma regression models
    !-----------------------------------------
    sll = 0.d0
    do t = (m +1),n
       sll = sll + d_gamma(y(t), 2, (/mu(t), nu/), .true.)
    end do
    return
  end function llk_gamma

  subroutine dllk_gamma(m,n,y,n1,mut,skipmu,n2,nut,skipnu,dllmu,dllnu,argsD)
    !-------------------------------------------------------------------------
    !
    !  Calculates dl/dmu and dl/dnu for gamma regression models
    !
    !   n1 = 1, if mu is fixed and n1 = n if mu is time varying
    !   n2 = 1, if mu is fixed and n2 = n if nu is time varying    
    !
    !-------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n, n1, n2, skipmu, skipnu
    real(dp), intent(in) :: y(n), mut(n1), nut(n2)
    real(dp), intent(out) :: dllmu(n*(1-skipmu)+skipmu)
    real(dp), intent(out) :: dllnu(n*(1-skipnu)+skipnu)
    type(argsDist), intent(in) :: argsD  ! dummy, for compatibility      
    real(dp) :: nu(n), mu(n), dgn(n)
    integer :: t
    real(dp) :: ym

    dllmu = argsD%arg1! Dummy to avoid compiler warning
    dllmu = 0.d0
    dllnu = 0.d0
    if(n1 + n2 == 0) return

    mu = mut(n1)
    nu = nut(n2)
    if(n1 > 1) mu = mut
    if(n2 > 1) nu = nut

    !-------------------------------------------------------
    ! dl/dnu
    !-------------------------------------------------------
    if(skipnu == 1) goto 200
    dgn = digamma(nu(n2))
    do t = (m +1), n
       ! updates only if n2 > 1
       if(n2 > 1) dgn(t) = digamma(nu(t))
       ym = y(t)/mu(t)
       dllnu(t) = 1.d0 - dgn(t) + log(nu(t)*ym) - ym
    end do

200 continue
    !-------------------------------------------------------------
    ! dl/dmu
    !--------------------------------------------------------------
    if(skipmu == 1) goto 300
    do t = (m +1), n
       dllmu(t) = nu(t)/mu(t) * (y(t)/mu(t)  - 1.d0)
    end do

300 continue
    return
  end subroutine dllk_gamma

  subroutine Ed2llk_gamma(m,n,n1,mut,i1,n2,nut,i2,E,argsD)
    !-------------------------------------------------------------------------
    !
    !  Calculates
    !           E(dl2/dmu2), E(dl2/dmudnu) and E(dl2/dnu2) 
    !  for gamma regression models
    !
    !   n1 = 1, if mu is fixed and n1 = n if mu is time varying
    !   n2 = 1, if mu is fixed and n2 = n if nu is time varying    
    !
    !-------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n, n1, n2, i1, i2
    real(dp), intent(in) :: mut(n1), nut(n2)
    real(dp), intent(out) :: E(n,max(1,i1+i2+i1*i2))
    type(argsDist), intent(in) :: argsD  ! dummy, for compatibility      
    real(dp) ::  nu(n), mu(n), psi
    integer :: t, ifault

    E = argsD%arg1! Dummy to avoid compiler warning
    E = 0.d0
    if(i1 + i2 == 0) return

    if(n1*n2 == 1) then
       !---------------------------------------
       !    mu and nu don't vary with time
       !---------------------------------------
       if(i1 == 1) E(:,i1) = nut(1)/mut(1)**2       
       if(i2 == 1) then
          if(i1 == 1) E(:,i1+i1*i2) = 0.d0
          E(:,i1+i2+i1*i2) = trigamma(nut(1), ifault) - 1/nut(1)
       end if
       return
    end if

    !---------------------------------------
    !    mu and/or nu vary with time
    !---------------------------------------
    mu = mut(n1)
    nu = nut(n2)
    if(n1 > 1) mu = mut
    if(n2 > 1) nu = nut

    psi = 0.d0
    if(n2*i2 == 1) psi =  trigamma(nu(1), ifault) 
    do t = (m+1),n       
       if(i1 == 1) E(t,i1) = nu(t)/mu(t)**2       
       if(i2 == 1) then
          if(i1 == 1) E(t,i1+i1*i2) = 0.d0
          if(n2 > 1) psi =  trigamma(nu(t), ifault)           
          E(t,i1+i2+i1*i2) = psi - 1/nu(t)
       end if
    end do

    return    
  end subroutine Ed2llk_gamma
  !--------------------------- GAMMA (START) ----------------------

  function llk_uw(m, n, y, mu, lambda, argsD) result(sll)
    !-----------------------------------------------------------------------
    ! Computes the sum(log(f(y| mu, nu))) for Unit-Weibull regression models
    ! argsD contains the quantile of interest rho, 
    ! for which mu is the rho-th quantile 
    !-----------------------------------------------------------------------
    implicit none    
    integer, intent(in) :: m, n
    real(dp), intent(in) :: y(n), mu(n), lambda
    type(argsDist), intent(in) :: argsD
    real(dp) :: sll
    integer :: t
    !-----------------------------------------------------
    ! log-likelihood for Unit-Weibull regression models
    !-----------------------------------------------------
    sll = 0.d0
    do t = (m+1),n
       sll = sll + duw(y(t), 3, (/mu(t), lambda, argsD%arg1/), .true.)
    end do
    return
  end function llk_uw


  subroutine dllk_uw(m,n,y,n1,mut,skipmu,n2,nut,skipnu,dllmu,dllnu, argsD)
    !-------------------------------------------------------------------------
    !
    !  Calculates dl/dmu and dl/dnu for Unit-Weibull regression models
    !
    !  Variable names are similar to the Kuma and Beta cases, for simplicity
    !  argsD contains rho, which is a constant in the model.
    !
    !-------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n, n1, n2, skipmu, skipnu
    real(dp), intent(in) :: y(n), mut(n1), nut(n2)
    real(dp), intent(out) :: dllmu(n*(1-skipmu)+skipmu)
    real(dp), intent(out) :: dllnu(n*(1-skipnu)+skipnu)
    type(argsDist), intent(in) :: argsD
    real(dp) :: At(n), num(n), denum(n), rho
    real(dp) :: mu(n), lambda(n)
    integer :: t

    dllmu = 0.d0
    dllnu = 0.d0
    if(n1 + n2 == 0) return

    mu = mut(n1)
    lambda = nut(n2)
    rho = argsD%arg1

    if(n1 > 1) mu = mut
    if(n2 > 1) lambda = nut

    !-------------------------------------------------------   
    ! At
    !-------------------------------------------------------
    do t = (m + 1), n
       At(t) = log(y(t))/log(mut(t))
    end do

    !-------------------------------------------------------
    ! dl/dnu
    !-------------------------------------------------------
    if(skipnu == 1) goto 200

    do t = (m + 1), n
       dllnu(t) = 1/lambda(t) + (1 + log(rho)*At(t)**lambda(t))*log(At(t))
    end do


200 continue
    !-------------------------------------------------------------
    ! dl/dmu
    !--------------------------------------------------------------
    if(skipmu == 1) goto 300
    do t = (m + 1), n
       num(t) = lambda(t)*(1+log(rho)*At(t)**lambda(t))
       denum(t) = mut(t)*log(mut(t))
       dllmu(t) = - num(t)/denum(t)
    end do

300 continue
    return
  end subroutine dllk_uw

  subroutine Ed2llk_uw(m,n,n1,mut,i1,n2,nut,i2,E, argsD)
    !-------------------------------------------------------------------------
    !
    !   Calculates
    !           E(dl2/dmu2), E(dl2/dmudnu) and E(dl2/dnu2) 
    !   for Unit-Weibull regression models. These are the columns 
    !   of E in the code
    !
    !   n1 = 1, if mu is fixed and n1 = n if mu is time varying
    !   n2 = 1, if nu is fixed and n2 = n if nu is time varying   
    !   i1/i2 = 0, if mu/nu are fixed values provided by the user
    !
    !-------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: m, n, n1, n2, i1, i2
    real(dp), intent(in) :: mut(n1), nut(n2)
    real(dp), intent(out) :: E(n,max(1,i1+i2+i1*i2))
    type(argsDist), intent(in) :: argsD    
    real(dp) :: cons0, cons1, cons2, nu(n), mu(n), lr
    integer :: t
    real(dp), parameter :: em = 0.57721566490153286061d0

    E = 0.d0
    if(i1 + i2 == 0) return

    mu = mut(n1)
    nu = nut(n2)
    if(n1 > 1) mu = mut
    if(n2 > 1) nu = nut

    lr = log(-log(argsD%arg1))
    cons0 = em + lr
    cons1 = cons0 - 1
    cons2 = lr*(lr + 2*em - 2) + pi**2/6 + (em - 2)*em

    if(n1*n2 == 1) then
       !---------------------------------------
       !    mu and lambda are not time varying
       !---------------------------------------
       if(i1 == 1) E(:,i1) = (nu(1)/(mu(1)*log(mu(1))))**2
       if(i1*i2 == 1) E(:,i1+i1*i2) = cons1/(mu(1)*log(mu(1)))
       if(i2 == 1) E(:,i1+i2+i1*i2) = (1.d0 + cons2)/nu(1)**2
    else
       !---------------------------------------------
       ! at least one varies with time
       !---------------------------------------------             
       do t = (m + 1), n
          if(i1 == 1) E(t,i1) = (nu(t)/(mu(t)*log(mu(t))))**2
          if(i1*i2 == 1) E(t,i1+i1*i2) = cons1/(mu(t)*log(mu(t)))
          if(i2 == 1) E(t,i1+i2+i1*i2) = (1.d0 + cons2)/nu(t)**2
       end do
    end if


    return    
  end subroutine Ed2llk_uw

  !--------------------------- Unit Weibull (END) --------------------

end module base


subroutine linkR(link, a, ylim, n, ilk, y, lk, gy, dl, dlink)
  !------------------------------------------------------------------------
  !  Interface for R
  !  Evaluates:
  !     - the link function at y, if lk = 1
  !     - the inverse of the link function at gy, if ilk = 1
  !     - the derivative of the link function at y, if dl = 1
  !
  !  Input:
  !     link = an integer indicating the type of link
  !            0 = linear      --> g(y) = a*y
  !            1 = logit       --> g(y) = log((y-ymin)/(ymax-y))
  !            2 = logarithm   --> g(y) = log(y)
  !            3 = loglog      --> g(y) = log(-log((y - ymin)/(ymax - ymin)))
  !            4 = c-loglog    --> g(y) = log(-log(1.d0 - (y - ymin)/(ymax - ymin)))
  !
  !     a = constant for linear link (will be ignored for other links)
  !
  !     ylim = vector of size 2 containing the lower and upper limits for y
  !
  !     n = sample size
  !
  !     ilk = 1 indicates that the iverse of the link must be calculated
  !           using gy
  !
  !     y = the values of the variable y
  !
  !     lk = 1 indicates that the link must be evaluated using y
  !
  !     gy = the values of g(y)
  !
  !     dl = 1 indicates that the derivative of the link must be evaluated
  !          using y
  !
  !     dlink = the values of dg(y)/dy
  !------------------------------------------------------------------------
  use base
  implicit none
  integer :: n, link, lk, dl, ilk
  real(dp) :: y(n), a, ylim(2)
  real(dp) :: gy(max(1,n*lk)), dlink(max(1,n*dl))
  type(argslink) :: args
  integer :: i
  args%link = link
  args%a = a    
  args%lower = ylim(1)
  args%upper = ylim(2)
  do i = 1,n
     if(lk == 1) gy(i) = linkfun(y(i), args)
     if(ilk == 1) y(i) = linkinv(gy(i), args)
     if(dl == 1) dlink(i) = diflink(y(i), args)
  end do
  return
end subroutine linkR

subroutine btsrpredictR(n, y, ylower, yupper, gy, nreg, xreg, escale, error, &
     nnew, xnew, ynew, link, npar, par, fixa, alpha, fixB, flagsb, fvbeta, &
     p, fixar, flagsar, fvar, xregar, q, fixma, flagsma, fvma, fixd, d, fixnu, nu, &
     inf) 
  !-----------------------------------------------------------------------------------------
  !
  !  Subrotuine used for prediction.
  !  The values of y, gy, eta, error, 
  !
  !-----------------------------------------------------------------------------------------
  use base
  implicit none
  integer:: n, nreg, link(2), inf, escale, nnew
  integer :: npar, fixa, fixB, fixar, fixma, fixd, fixnu
  integer :: p, q, xregar
  integer :: flagsb(max(1,fixB)),flagsar(max(1,fixar)), flagsma(max(1,fixma))
  real(dp) :: y(n), xreg(n, max(1, nreg)), ylower, yupper
  real(dp) :: par(npar), alpha, d, nu, fvbeta(max(1,fixB))
  real(dp):: fvar(max(1,fixar)), fvma(max(1,fixma))    
  real(dp) :: gy(n), error(n), xnew(nnew, max(1,nreg)), ynew(nnew)
  real(dp) :: xstart(max(1, nreg))
  type(argsModel) :: model
  real(dp), allocatable :: vc(:)

  !---------------------------------------------
  !  If d = 0 uses inf = q
  !---------------------------------------------
  inf = max(inf, q)
  if(d == 0.0d0) inf = q

  xstart = 0.d0  ! dummy

  !----------------------------------------------------------------
  ! setting the link parameters
  !----------------------------------------------------------------
  if(.not. allocated(model%argsL)) allocate(model%argsL(2))        
  model%argsL(1:2)%link = link
  model%argsL(1:2)%lower = ylower
  model%argsL(1:2)%upper = yupper         

  !-------------------------------------------------------------
  ! allocating the time series and parameters
  ! setting fixed values/lags of parameters
  !-------------------------------------------------------------
  call allocate_model(model, n, y, gy, nreg, xreg, xstart,1-fixnu, nu, 1-fixa, alpha, &
       nreg-fixb, flagsb, fvbeta, xregar, p, p-fixar, flagsar, fvar, q, q-fixma, flagsma, fvma,&
       1-fixd, d, inf, 0)

  !-------------------------------------------------------------
  ! setting the values of eta and error (from 1 to n)
  !-------------------------------------------------------------
  model%error = error
  model%cts(1)%eta = 0.d0
  model%escale = escale

  !-----------------------------------------------------
  !  Setting the parameter's non-fixed values:
  !  alpha, beta, ar, ma and d
  !-----------------------------------------------------
  call safe_allocate(vc, 0, model%inf(1))
  call start_par(par, model, vc, 1)

  !-----------------------------------------------------
  !  predicted values
  !-----------------------------------------------------      
  call mu_forecast(model, vc, nnew, xnew, ynew)

  return    
end subroutine btsrpredictR
