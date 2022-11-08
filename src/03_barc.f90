module barc
  use main_mod
  use base
  implicit none

contains

  function map_T(x, r, theta, mtype) result(Tx)
    implicit none
    real(dp), intent(in) :: x
    integer, intent(in) :: r, mtype
    real(dp), intent(in) :: theta(r)
    real(dp) :: Tx

    Tx = 0.d0
    select case(mtype)
    case (1)
       ! (kx)(mod 1), k must be an integer greater or equal than 2
       Tx = theta(1)*x
       Tx = Tx - int(Tx)
    case (2)
       ! Rafael's map. 0 <= theta <= 1
       if(x < theta(1)) then
          Tx = x/theta(1)
       else
          Tx = theta(1)*(x - theta(1))/(1.d0 - theta(1))
       end if
    case (3)
       ! logistic map. 0 <= theta  <= 4
       Tx = theta(1)*x*(1.d0 - x)
    case (4)
       ! Manneville-Pomeau. 0 < theta < 1
       Tx = x + x**(1+theta(1))
       Tx = Tx - int(Tx)
    case (5)
       ! Lasota-Mackey's map. No theta
       if(x <= 0.5d0) then
          Tx = x/(1.d0 - x)
       else
          Tx = 2.d0*x - 1.d0
       end if
    end select

  end function map_T

  subroutine start_par_barc(par, model)
    !--------------------------------------------------------
    !
    !             Parameter Initialization
    !
    ! Uses the subroutine from base module to allocate
    ! alpha, beta and phi. Allocates thetaT, u0 and nu
    !
    !---------------------------------------------------------
    implicit none
    type(argsModel) :: model
    real(dp) :: par(model%npar(1))
    integer :: n1, n2

    call start_par(par, model, 1)

    ! fixing nu and setting the initial values for u0 (if needed) and theta_T
    n1 = model%alpha(1)%fit + model%beta(1)%fit + model%ar(1)%fit
    n2 = n1

    ! ThetaT
    n1 = n2 + 1
    n2 = n2 + model%thetaT%fit
    if(n2 >= n1) model%thetaT%par = par(n1:n2)

    ! fixing nu
    n1 = n2 + 1
    n2 = n2 + model%nu%fit
    if(n2 >= n1) model%nu%par = par(n1:n2)

    ! u0 is the last value in the parameter's vector
    if(model%u0%fit == 1) model%u0%par = par(model%npar(1))

    return    
  end subroutine start_par_barc

  subroutine mu_calc_barc(model)
    !--------------------------------------------------------
    !
    !            Recurrence for mu(t)
    !
    ! Uses the subroutine from base module to calculate
    ! the recurrence assuming a MA(0) and d = 0. Here the 
    ! chaotic part is added
    !
    !---------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model    
    real(dp) :: vc(1)
    integer :: t

    vc = 0.d0
    ! error scale makes no difference here:
    call mu_calc(model%n,model%y,model%gy, model%ystart,&
         model%cts(1)%nreg,model%cts(1)%xreg, model%cts(1)%xstart,&         
         model%cts(1)%ut,model%cts(1)%eta,model%error,0,&         
         model%alpha(1)%par(1),model%beta(1)%par,model%ar(1)%length, model%ar(1)%par, &
         model%cts(1)%xregar,0,vc,model%m, model%argsL(1:2))

    !-----------------------------------------------
    ! calculating T**{t-1}(u0) and adding to eta_t
    !-----------------------------------------------
    model%cts(1)%orbit(1) = model%u0%par(1)
    model%cts(1)%eta(1) = model%cts(1)%eta(1) + linkfun(model%cts(1)%orbit(1),model%argsL(3))
    model%cts(1)%ut(1) = linkinv(model%cts(1)%eta(1),model%argsL(1))
    ! to avoid the bounds:
    if(model%cts(1)%ut(1) <= 0.d0) then
       model%cts(1)%ut(1) = epsilon(1.d0)
       model%cts(1)%eta(1) = linkfun(model%cts(1)%ut(1), model%argsL(1))
    elseif(model%cts(1)%ut(1) >= 1.d0) then
       model%cts(1)%ut(1) = 1.d0 - epsilon(1.d0)
       model%cts(1)%eta(1) = linkfun(model%cts(1)%ut(1), model%argsL(1))
    end if

    do t = 2, model%n
       ! T**{t-1}(u0)
       model%cts(1)%orbit(t) = map_T(model%cts(1)%orbit(t-1), model%ThetaT%length, model%ThetaT%par, model%map)
       ! eta + h(T**{t-1}(u0))
       model%cts(1)%eta(t) = model%cts(1)%eta(t) + linkfun(model%cts(1)%orbit(t),model%argsL(3))
       ! mut = g**{-1}(eta_t)
       model%cts(1)%ut(t) = linkinv(model%cts(1)%eta(t), model%argsL(1))
       ! to avoid the bounds:
       if(model%cts(1)%ut(t) <= 0.d0) then
          model%cts(1)%ut(t) = epsilon(1.d0)
          model%cts(1)%eta(t) = linkfun(model%cts(1)%ut(t), model%argsL(1))
       elseif(model%cts(1)%ut(t) >= 1.d0) then
          model%cts(1)%ut(t) = 1.d0 - epsilon(1.d0)
          model%cts(1)%eta(t) = linkfun(model%cts(1)%ut(t), model%argsL(1))
       end if
    end do

    if(model%escale == 0) then
       model%error = model%y - model%cts(1)%ut
    else
       do t = 1,model%n             
          if(model%argsL(1)%link /= model%argsL(2)%link) then       
             model%error(t) = linkfun(model%y(t), model%argsL(1)) - model%cts(1)%eta(t)
          else
             model%error(t) = model%gy(t) - model%cts(1)%eta(t)       
          end if
       end do
    end if

    return
  end subroutine mu_calc_barc

  subroutine mu_forecast_barc(model, nnew, xhat, yhat, That)
    !--------------------------------------------------------
    !
    !   Forecast for a BARC model
    !
    !---------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    integer, intent(in) :: nnew
    real(dp), intent(in) :: xhat(nnew, max(1, model%beta(1)%length))
    real(dp), intent(out) :: yhat(nnew), That(nnew)
    real(dp) :: Ts(0:nnew)
    real(dp) :: XB((model%n+1-model%ar(1)%length):(model%n+nnew))
    real(dp) :: gyhat(nnew)
    real(dp) :: gy((model%n+1-model%ar(1)%length):(model%n+nnew))
    integer :: t, i

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

    Ts(0) = model%cts(1)%orbit(model%n)
    do t = 1,nnew
       ! ghat(y) = a + x*b + h(T^{t-1}(u0))
       Ts(t) = map_T(Ts(t-1), model%ThetaT%length, model%ThetaT%par, model%map)
       gyhat(t) =  model%alpha(1)%par(1) + XB(model%n+t) + linkfun(Ts(t),model%argsL(3))
       ! ghat(y) = a + x*b + h(T^{t-1}(u0)) + AR
       do i = 1,model%ar(1)%length   ! if p < 1 ignores the loop
          gyhat(t) = gyhat(t) + model%ar(1)%par(i)*gy(model%n+t-i)
          if(model%cts(1)%xregar == 1) gyhat(t) = gyhat(t) - model%ar(1)%par(i)*XB(model%n+t-i)
       end do
       ! yhat = g1^{-1}(gyhat)
       yhat(t) = linkinv(gyhat(t), model%argsL(1))
       ! g2(y) to be used in the AR recursion
       if(model%argsL(1)%link == model%argsL(2)%link) then
          gy(model%n+t) = gyhat(t)
       else
          gy(model%n+t) = linkinv(yhat(t),model%argsL(2))
       end if
    end do
    That = Ts(1:nnew)
    return
  end subroutine mu_forecast_barc

  subroutine get_model_barc(model, n, y, gy, ystart, nreg, xreg, xregar, xstart, escale,&
       link, ah, npar, par, fixa, alpha, fixB, flagsb, fvbeta,&
       p, fixar, flagsar, fvar, r, fixt, flagst, fvt, fixnu, nu, fixu0, u0,  &
       llk, sco, info, map)
    !-----------------------------------------------------------------------------------------
    !
    !  Subrotuine used to pass the values enter in the main program to the user defined
    !  variable that will be passed to the  generic subroutines
    !
    !-----------------------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    integer, intent(in) :: n, nreg, link(3), llk, sco, info, map, xregar, escale
    integer, intent(in) :: npar, fixa, fixB, fixar, fixt, fixnu, p, r, fixu0
    integer, intent(in) :: flagsb(max(1,fixB)),flagsar(max(1,fixar)), flagst(max(1,fixt))
    real(dp), intent(in) :: y(n), xreg(n, max(1, nreg)), u0, ah
    real(dp), intent(in) :: par(npar), alpha, nu, fvbeta(max(1,fixB))
    real(dp), intent(in):: fvar(max(1,fixar)), fvt(max(1,fixt))
    real(dp), intent(inout) :: gy(n), ystart, xstart(max(1,nreg))

    call get_model(model, n, y, 0.d0, 1.d0, ystart,  gy, nreg, xreg, xstart, link(1:2),&
         escale, 1, npar, par, fixa, alpha, fixB, flagsb, fvbeta,&
         p, fixar, flagsar, fvar,  xregar, 0, 0, (/0/), (/0.d0/), 1, 0.d0, &
         fixnu, nu, 0, llk, sco, info, 0)

    ! chaotic map
    model%map = map

    !----------------------------------------------------------------
    ! setting the link parameters
    ! setting the link parameters for h(T**t(u0))
    !----------------------------------------------------------------
    if(allocated(model%argsL)) deallocate(model%argsL)
    allocate(model%argsL(3))
    model%argsL(:)%link = link
    model%argsL(:)%lower = 0.d0
    model%argsL(:)%upper = 1.d0
    model%argsL(3)%a = ah

    !-------------------------------------------------------------
    ! allocating the orbit
    !-------------------------------------------------------------
    call safe_allocate(model%cts(1)%orbit, n)

    !-------------------------------------------------------------
    ! theta_T and u0
    !-------------------------------------------------------------
    call allocate_parvec(model%thetaT, r, fixt, flagst, fvT)
    call allocate_parvec(model%u0, 1, fixu0, (/1/), (/u0/))
    model%npar(1) = model%npar(1) + (r - fixT) + (1-fixu0)

    !------------------------------------------------------
    !  Setting the parameter's fixed values
    !-----------------------------------------------------
    call start_par_barc(par, model)

    if(sco + info == 0 ) return
    !-------------------------------------------------------------------------------------------
    ! allocating score-vector, Information matrix and related matrices
    !-------------------------------------------------------------------------------------------
    !!  PRECISA ARRUMAR AQUI QUANDO FIZER AS DERIVADAS TEORICAS 
    !!  ALOCAR O QUE NÃO FOI ALOCADO NA FUNÇÃO BASE           

    return    
  end subroutine get_model_barc


  subroutine return_model_barc(model, mu, eta, error, Ts)
    implicit none
    type(argsModel), intent(in) :: model
    real(dp), intent(out), dimension(model%n) :: mu, eta, error, Ts

    mu = model%cts(1)%ut
    eta = model%cts(1)%eta
    error = model%error(1:model%n)
    Ts = model%cts(1)%orbit
    return
  end subroutine return_model_barc

  subroutine U_barc_numeric(model, npar, par, U)
    implicit none
    integer, intent(in) :: npar
    type(argsModel), intent(inout) :: model
    real(dp), intent(in) :: par(npar)
    real(dp), intent(out) :: U(npar)
    integer :: i
    real(dp), parameter :: eps = 1.d-4
    real(dp) :: par1(npar), par2(npar), f1, f2

    !-----------------------------------------
    ! Score vector - numerical derivative
    !-----------------------------------------
    U(:) = 0.d0
    do i = 1, npar
       par1 = par
       par2 = par
       par1(i) = par1(i)+eps
       par2(i) = par2(i)-eps
       call start_par_barc(par1, model)    
       call mu_calc_barc(model)
       f1 = llk_beta(model%m, model%n, model%y, model%cts(1)%ut, model%nu%par(1), model%argsD)
       call start_par_barc(par2, model)    
       call mu_calc_barc(model)
       f2 = llk_beta(model%m, model%n, model%y, model%cts(1)%ut, model%nu%par(1), model%argsD)
       U(i) = (f1 - f2)/(2.*eps)          
    end do

    ! Restoring parameters
    call start_par_barc(par, model)    
    return
  end subroutine U_barc_numeric

  subroutine K_barc_numeric(model, npar, par, K)
    implicit none
    integer, intent(in) :: npar
    type(argsModel), intent(inout) :: model
    real(dp), intent(in)  :: par(npar)
    real(dp), intent(out) :: K(npar, npar)
    real(dp), parameter :: eps = 1.d-4
    real(dp) :: np1(npar), np2(npar), np3(npar), np4(npar), h1, h2, h3, h4
    integer :: i, j

    !-------------------------
    ! Information matrix
    !-------------------------    
    do i = 1,npar
       do j = 1,i
          np1 = par; np2 = par; np3 = par; np4 = par
          np1(i) = np1(i) + eps; np1(j) = np1(j) + eps
          np2(i) = np2(i) + eps; np2(j) = np2(j) - eps
          np3(i) = np3(i) - eps; np3(j) = np3(j) + eps
          np4(i) = np4(i) - eps; np4(j) = np4(j) - eps        
          call start_par_barc(np1, model)    
          call mu_calc_barc(model)
          h1 = llk_beta(model%m, model%n, model%y, model%cts(1)%ut, model%nu%par(1), model%argsD)
          call start_par_barc(np2, model)    
          call mu_calc_barc(model)
          h2 = llk_beta(model%m, model%n, model%y, model%cts(1)%ut, model%nu%par(1), model%argsD)
          call start_par_barc(np3, model)    
          call mu_calc_barc(model)
          h3 = llk_beta(model%m, model%n, model%y, model%cts(1)%ut, model%nu%par(1), model%argsD)
          call start_par_barc(np4, model)    
          call mu_calc_barc(model)
          h4 = llk_beta(model%m, model%n, model%y, model%cts(1)%ut, model%nu%par(1), model%argsD)
          K(i, j) = (h1 -h2-h3 +h4 ) / (4.*eps*eps)
          K(j,i) = K(i, j)                             
       end do
    end do
    ! Information is -Hess.
    K = -K

    ! Restoring parameters
    call start_par_barc(par, model)   
    return
  end subroutine K_barc_numeric


  subroutine loglik_barc(model, npar, par, sll, U)
    !------------------------------------------------------------------
    !
    !   Log-likelihood: BARC model
    !
    !------------------------------------------------------------------
    implicit none
    integer, intent(in) :: npar
    real(dp), intent(in) :: par(npar)
    type(argsModel), intent(inout) :: model
    real(dp), intent(out) :: sll, U(npar)

    ! Initializing parameters 
    call start_par_barc(par, model)    

    U = 0.d0
    if(model%sco == 0) goto 100

    ! ordem trocada por causa da derivada numérica
    !-----------------------------------------
    !    Score vector
    !-----------------------------------------    
    call U_barc_numeric(model, npar, par, U)
    U = -U   

100 continue   

    ! Calculating recursions that do not depend on the distribution
    call mu_calc_barc(model)

    !-----------------------------------------
    ! log-likelihood for BARC model
    !-----------------------------------------
    sll = llk_beta(model%m, model%n, model%y, model%cts(1)%ut, model%nu%par(1), model%argsD)
    sll = -sll
    !if(sll < -Huge(1.d0)) sll = -Huge(1.d0)
    !if(sll > Huge(1.d0)) sll = Huge(1.d0)
    return
  end subroutine loglik_barc


  subroutine loglik_barc_nelder(model, npar, par, sll)
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

    model%llk = 1
    model%sco = 0
    call loglik_barc(model, npar, par_aux, sll, U)

    return
  end subroutine loglik_barc_nelder

  subroutine loglik_barc_lbfgsb(model, npar, par, sll, U)
    !------------------------------------------------------------------
    !
    !   Subroutine to be used in L-BFGS-B optimization subroutine
    !
    !------------------------------------------------------------------
    implicit none
    integer, intent(in) :: npar
    real(dp), intent(in) :: par(npar)
    type(argsModel), intent(inout) :: model
    real(dp), intent(out) :: sll, U(npar)

    model%llk = 1
    model%sco = 1
    call loglik_barc(model, npar, par, sll, U)

    return
  end subroutine loglik_barc_lbfgsb

end module barc

subroutine barcR(n, y, gy, ystart, nreg, xreg, xstart, mu, eta, error, escale, Ts, nnew, xnew, ynew, Tnew, &
     link, ah, map, npar, par, fixa, alpha, fixB, flagsb, beta, p, fixar, flagsar, phi, xregar,&
     r, fixt, flagst, thetaT, fixnu, nu, fixu0, u0, llk, sll, sco, U, info, K)
  use barc
  implicit none
  integer :: n, nreg, link(3), map, nnew, xregar, escale
  real(dp) :: y(n), gy(n), xreg(n, max(1, nreg)), Ts(n), ah
  real(dp) :: mu(n), eta(n), error(n), ynew(max(1,nnew)), ystart
  real(dp) :: Tnew(max(1,nnew)), xnew(max(1,nnew),max(1, nreg))
  integer :: npar, fixa, fixB, fixar, fixt, fixnu, fixu0, p, r
  integer :: flagsb(max(1,fixB)), flagsar(max(1,fixar)), flagst(max(1,fixt))
  real(dp) :: par(npar), alpha, nu, beta(max(1,fixB))
  real(dp) :: phi(max(1,fixar)), thetaT(max(1,fixt)), u0
  integer :: llk, sco, info
  real(dp) :: sll,  U(max(1, npar*sco)), xstart(max(1,nreg))
  real(dp) :: K(max(1, npar*info), max(1, npar*info)), Utemp(npar)
  type(argsModel) :: model

  !---------------------------------------------------------------------------------------------------
  ! allocating matrices and vectors and setting variables fixed values
  !---------------------------------------------------------------------------------------------------
  call get_model_barc(model, n, y, gy, ystart, nreg, xreg, xregar, xstart, escale, link, ah,&
       npar, par, fixa, alpha, fixB, flagsb, beta, p, fixar, flagsar, phi, r, fixt, flagst, thetaT, &
       fixnu, nu, fixu0, u0, llk, sco, info, map)    

  ! ordem trocada por causa da hessiana numérica
  K = 0.d0
  if(info == 0)  goto 100

  call K_barc_numeric(model, npar, par, K)

100 continue          

  sll = 0.d0
  U = 0.d0
  ! calculates:
  !    mu, eta, Ts
  !    -sll: log-likelihood (if llk = 1)
  !    -U: score vector (if sco = 1)
  call loglik_barc(model, npar, par, sll, Utemp)

  call return_model_barc(model, mu, eta, error, Ts)
  sll = -sll
  if(sco == 1) U = -Utemp

  if(nnew == 0) return
  call mu_forecast_barc(model, nnew, xnew, ynew, Tnew)

  return      
end subroutine barcR


subroutine optimnelderbarcR(npar, par, nbd, lower, upper, n, y, gy, ystart, nreg, xreg, xstart, &
     mu, eta, error, escale, Ts, nnew, xnew, ynew, Tnew, link, ah, map, fixa, alpha, fixB, flagsb, beta,&
     p, fixar, flagsar, phi,  xregar, r, fixt, flagst, thetaT, fixnu, nu, fixu0, u0, sll, sco, U, info, K,&
     iprint, stopcr, maxit, neval, conv)
  use barc
  implicit none
  !--------------------------------------------
  ! iprint < 0 = no print
  ! stopcr = stopping critereon  (1.d-4)
  !---------------------------------------------
  integer :: n, nreg, link(3), map, nnew, escale
  integer :: npar, fixa, fixB, fixar, xregar, fixt, fixnu, fixu0, p, r
  integer :: flagsb(max(1,fixB)), flagsar(max(1,fixar)), flagst(max(1,fixt))
  integer ::  nbd(npar), neval, sco, info, maxit, iprint, conv    
  real(dp) :: y(n), gy(n), xreg(n, max(1, nreg)), Ts(n), ah
  real(dp) :: mu(n), eta(n), error(n), ynew(max(1,nnew)), ystart
  real(dp) :: Tnew(max(1,nnew)), xnew(max(1,nnew),max(1, nreg))
  real(dp) :: par(npar), alpha, nu, beta(max(1,fixB)), xstart(max(1,nreg))
  real(dp) :: phi(max(1,fixar)), thetaT(max(1,fixt)), u0
  real(dp) :: U(max(1, npar*sco)), K(max(1, npar*info),max(1, npar*info))
  real(dp) :: sll, lower(npar), upper(npar),stopcr
  type(argsModel) :: model
  real(dp) :: Utemp(npar)

  !---------------------------------------------------------------------------------------------------
  ! allocating matrices and vectors and setting variables fixed values
  !---------------------------------------------------------------------------------------------------
  call get_model_barc(model, n, y, gy, ystart, nreg, xreg, xregar,xstart, escale, link, ah, &
       npar, par, fixa, alpha, fixB, flagsb, beta, p, fixar, flagsar, phi, &
       r, fixt, flagst, thetaT, fixnu, nu, fixu0, u0, 1, 0, 0, map)

  !-------------------------------------------------------------------------------
  ! Nelder-Mead
  !-------------------------------------------------------------------------------
  call optim(loglik_barc_nelder, model, npar, par, nbd, lower, upper, sll, iprint,&
       stopcr, maxit, neval, conv)

  K = 0.d0
  if(info == 0)  goto 200
  call K_barc_numeric(model, npar, par, K)

200 continue      
  !---------------------------------------------------------------------------------------------------
  ! Reports the final model
  !---------------------------------------------------------------------------------------------------
  model%sco = sco   
  call loglik_barc(model, npar, par, sll, Utemp)

  call return_model_barc(model, mu, eta, error, Ts)
  sll = -sll
  U = 0.d0
  if(sco == 1) U = -Utemp

  if(nnew == 0) return
  call mu_forecast_barc(model, nnew, xnew, ynew, Tnew)

  return      
end subroutine optimnelderbarcR


subroutine optimlbfgsbbarcR(npar, par, nbd, lower, upper, n, y, gy, ystart, nreg, xreg,xstart,&
     mu, eta, error, escale, Ts, nnew, xnew, ynew, Tnew, link, ah, map, fixa, alpha, fixB, flagsb, beta,&
     p, fixar, flagsar, phi,  xregar, r, fixt, flagst, thetaT,  fixnu, nu,  fixu0, u0, sll, U, info, K, &
     iprint, factr, pgtol, maxit, neval, conv) 
  use barc
  implicit none
  integer :: n, nreg, link(3), map, conv, nnew, escale
  integer :: npar, fixa, fixB, fixar, xregar, fixt, fixnu, p, r, fixu0 
  integer :: flagsb(max(1,fixB)), flagsar(max(1,fixar)), flagst(max(1,fixt))
  integer ::  nbd(npar), neval, info,  maxit, iprint
  real(dp) :: y(n), gy(n), xreg(n, max(1, nreg)), ah, ystart
  real(dp) :: mu(n), eta(n), error(n), Ts(n), u0, ynew(max(1,nnew))
  real(dp) :: Tnew(max(1,nnew)), xnew(max(1,nnew),max(1, nreg))
  real(dp) :: par(npar), alpha,  nu, beta(max(1,fixB)), xstart(max(1,nreg))
  real(dp) :: phi(max(1,fixar)), thetaT(max(1,fixt))
  real(dp) :: sll, U(npar), K(max(1, npar*info),max(1, npar*info))
  real(dp) :: factr, pgtol, lower(npar), upper(npar)
  character(len = 60) :: conve
  type(argsModel) :: model

  !---------------------------------------------------------------------------------------------------
  ! allocating matrices and vectors and setting variables fixed values
  !---------------------------------------------------------------------------------------------------
  call get_model_barc(model, n, y, gy, ystart, nreg, xreg, xregar, xstart, escale, link, ah, &
       npar, par, fixa, alpha, fixB, flagsb, beta,  p, fixar, flagsar, phi, r, fixt, flagst, thetaT,&
       fixnu, nu, fixu0, u0, 1, 0, 0, map)

  !---------------------------------------------------------------------------------------------------
  ! L-BFGS-B
  !---------------------------------------------------------------------------------------------------
  call optim(loglik_barc_lbfgsb, model, npar, par, nbd, lower, upper,&
       sll, U, iprint, factr, pgtol, maxit, neval, conv, conve)

  K = 0.d0
  if(info == 0)  goto 200
  call K_barc_numeric(model, npar, par, K)

200 continue
  !---------------------------------------------------------------------------------------------------
  ! Reports the final model
  !---------------------------------------------------------------------------------------------------
  model%llk = 1
  model%sco = 1
  call loglik_barc(model, npar, par, sll, U)    

  sll = -sll
  U = -U  
  call return_model_barc(model, mu, eta, error, Ts)

  if(nnew == 0) return
  call mu_forecast_barc(model, nnew, xnew, ynew, Tnew)
  return      
end subroutine optimlbfgsbbarcR


subroutine simbarcR(n, burn, nu, alpha, nreg, beta, p, phi, r, theta, u0, map, link, ah, &
     xreg, xregar, yt, ystart, xstart, mu, eta,  error, escale, Ts, ns, seed, rngtype, rev) 
  use barc
  !#---------------------------------------------------------------------
  !#
  !#  Simulating a BARC model
  !#
  !#---------------------------------------------------------------------
  implicit none
  integer :: n, burn, nreg, p, r, map, link(3), ns
  integer :: seed(ns), rngtype, rev, xregar, escale
  real(dp) :: alpha, beta(max(1,nreg)), ah, nu
  real(dp) :: xreg(n+burn,max(1, nreg)), ystart, xstart(max(1,nreg))
  real(dp) :: phi(max(1,p)), u0, theta(max(1,r)), gy(n+burn)
  real(dp) :: yt(n+burn), mu(n+burn), eta(n+burn)
  real(dp) :: Ts(n+burn), error(n+burn)
  integer ::  t, i
  type(argslink) :: argsL(3)
  type(rng_t) :: rng
  real(dp) :: ytemp, xb, gytemp

  ! revision required
  rev = 1

  !# check for regressors
  if(nreg == 0) then
     xreg = 0.d0
     beta = 0.d0
  end if

  !# check for AR component
  if(p == 0) phi = 0.d0

  !# initializing variables in the model
  eta = 0.d0
  mu = 0.d0
  Ts = 0.d0
  yt = 0.d0
  gy = 0.d0
  error = 0.d0

  argsL(:)%link = link  ! g1(mu), g2(y), h(T(u))  
  argsL(3)%a = ah
  argsL(:)%lower = 0.d0
  argsL(:)%upper = 1.d0

  call rng_seed(rng, ns, seed, rngtype)

  ! starting values for g(yt)
  ytemp = 0.d0
  xb = 0.d0
  if(p > 1) then
     if(0.d0 < ystart .and. ystart < 1.d0) ytemp = linkfun(ystart, argsL(2))
     if(xregar == 1) xb = sum(xstart*beta)
  end if

  !# t = 1
  Ts(1) = u0
  eta(1) = alpha + sum(xreg(1,:)*beta) + (ytemp - xb)*sum(phi) + linkfun(Ts(1), argsL(3))

  mu(1) = linkinv(eta(1), argsL(1))
  if(mu(1) < 0.d0 .or. mu(1) > 1.d0) return
  yt(1) = rbeta(2, (/mu(1), nu/), rng)
  gy(1) = linkfun(yt(1), argsL(2))

  !# t > 1
  do t = 2, (n+burn)     
     !#-----------------
     !# T^{t-1}(u0)
     !#-----------------       
     Ts(t) =  map_T(Ts(t-1), r, theta, map)

     !#------------------------------------
     !# eta(t) = alpha + x*b + AR + h(T^{t-1}(u0))
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
     eta(t)  = eta(t) +  linkfun(Ts(t), argsL(3))

     !#---------------------------
     !# mu(t) = g^{-1}(eta(t))
     !#---------------------------
     mu(t) = linkinv(eta(t), argsL(1))
     if(mu(t) < 0.d0 .or. mu(t) > 1.d0) return

     if(mu(t) == 0.d0) then
        mu(t) = mu(t) + epsilon(1.d0)
        eta(t) = linkfun(mu(t), argsL(1))               
     elseif(mu(t) == 1.d0) then
        mu(t) = mu(t) - epsilon(1.d0)
        eta(t) = linkfun(mu(t), argsL(1))
     end if

     !#-----------------
     !# y(t) ~ beta
     !#-----------------
     yt(t) = rbeta(2, (/mu(t), nu/), rng)
     if(yt(t) == 0.d0) then
        yt(t) = yt(t) + epsilon(1.d0)
     elseif(yt(t) == 1.d0) then
        yt(t) = yt(t) - epsilon(1.d0)
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
end subroutine simbarcR

subroutine predictbarcR(n, y, gy, nreg, xreg, escale, error, Ts, nnew, xnew, ynew, Tnew,  &
     link, ah, map, npar, par, fixa, alpha, fixB, flagsb, fvbeta, p, fixar, flagsar, fvar, &
     xregar, r, fixt, flagst, fvT, fixnu, nu, fixu0, u0)
  use barc  
  !-----------------------------------------------------------------------------------------
  !
  !  Subrotuine used for prediction.
  !  The values of y, gy, eta, error, 
  !
  !-----------------------------------------------------------------------------------------
  implicit none
  integer:: n, nreg, link(3), escale, nnew
  integer :: npar, fixa, fixB, fixar, fixt,fixnu, fixu0
  integer :: p, r, xregar, map
  integer :: flagsb(max(1,fixB)),flagsar(max(1,fixar)), flagst(max(1,fixt))
  real(dp) :: y(n), xreg(n, max(1, nreg)), Ts(n), ah
  real(dp) :: par(npar), alpha, nu, fvbeta(max(1,fixB)), u0
  real(dp):: fvar(max(1,fixar)), fvT(max(1,fixt))    
  real(dp) :: gy(n), error(n), xnew(nnew, max(1,nreg)), ynew(nnew), Tnew(nnew)
  real(dp) :: xstart(max(1, nreg))
  type(argsModel) :: model

  xstart = 0.d0  ! dummy

  !----------------------------------------------------------------
  ! setting the link parameters
  ! setting the link parameters for h(T**t(u0))
  !----------------------------------------------------------------
  if(allocated(model%argsL)) deallocate(model%argsL)
  allocate(model%argsL(3))
  model%argsL(:)%link = link
  model%argsL(:)%lower = 0.d0
  model%argsL(:)%upper = 1.d0
  model%argsL(3)%a = ah

  ! chaotic map
  model%map = map

  !-------------------------------------------------------------
  ! allocating the time series and parameters
  ! setting fixed values/lags of parameters
  !-------------------------------------------------------------
  call allocate_model(model, n, y, gy, nreg, xreg, xstart,1-fixnu, nu, 1-fixa, alpha, &
       nreg-fixb, flagsb, fvbeta, xregar, p, p-fixar, flagsar, fvar, 0, 0, (/0/), (/0.d0/),&
       1, 0.d0, 0, 0)

  !-------------------------------------------------------------
  ! allocating the orbit
  !-------------------------------------------------------------
  call safe_allocate(model%cts(1)%orbit, n)

  !-------------------------------------------------------------
  ! theta_T and u0
  !-------------------------------------------------------------
  call allocate_parvec(model%thetaT, r, fixt, flagst, fvT)
  call allocate_parvec(model%u0, 1, fixu0, (/1/), (/u0/))
  model%npar(1) = model%npar(1) + (r - fixT) + (1-fixu0)

  !-------------------------------------------------------------
  ! setting the values of eta, error and orbit (from 1 to n)
  !-------------------------------------------------------------
  model%error = error
  model%cts(1)%eta = 0.d0
  model%escale = escale
  model%cts(1)%orbit = Ts      

  !-----------------------------------------------------
  !  Setting the parameter's non-fixed values:
  !  alpha, beta, ar, ma and d
  !-----------------------------------------------------
  call start_par_barc(par, model)

  !-----------------------------------------------------
  !  predicted values
  !-----------------------------------------------------      
  call mu_forecast_barc(model, nnew, xnew, ynew, Tnew)      

  return    
end subroutine predictbarcR

