module alloc
  use main_mod   ! user derived types
  !*******************************************************************************************
  !
  ! Allocation related subrotuines
  !
  !*******************************************************************************************
  ! February 2025
  !   - Added this module and moved some subroutines from other modules here
  !*******************************************************************************************
  implicit none


contains

  subroutine allocate_parvec(vec, length, fix, flags, fval)
    !*****************************************************************************************
    ! Subroutine used to allocate polynomials
    !*****************************************************************************************
    ! Input
    !  length: size of the vector
    !  fix   : number of fixed coefficients
    !  flags : the lags of the fixed coefficients
    !  fval  : values of the fixed parameters
    !
    ! Input / Output
    !  vec: vector of coefficients
    !*****************************************************************************************
    ! Last revision: July, 2023
    implicit none
    integer, intent(in) :: fix, length
    integer, intent(in) :: flags(max(1, fix))
    real(dp), intent(in) :: fval(max(1, fix))
    type(vec_parameter), intent(inout) :: vec
    integer :: dummy(max(1, length)), i

    vec%length = length
    vec%fit =  length - fix
    if (length == 0) return

    call safe_allocate(vec%cf, length)
    vec%cf = 0.d0
    dummy = 1

    !--------------------------------------------------------------
    ! fixed lags and fixed values
    !--------------------------------------------------------------
    if (fix > 0) then
       dummy(flags) = 0
       call safe_allocate(vec%flags, fix)
       vec%flags = flags
       vec%cf(vec%flags) = fval
    end if

    !--------------------------------------------------------------
    ! non-fixed lags
    !--------------------------------------------------------------
    if (length - fix > 0) then
       call safe_allocate(vec%lags, length - fix)
       vec%lags = pack([(i, i = 1, length)], dummy == 1)
    end if
    return
  end subroutine allocate_parvec

  subroutine allocate_conditional_ts(cts, n, nreg, xreg, xstart, part)
    !*****************************************************************************************
    ! Allocates the conditional variables and regressors for each part
    !*****************************************************************************************
    ! Input
    !  n      : sample size
    !  nreg   : number of regressors
    !  xreg   : matrix of regressors
    !  xstart : starting value for the regressors (when t < 1)
    !  part   : indicates which part of the model is being allocated (1 or 2)
    !
    ! Output
    !  cts: and object with the allocated vectors and matrices
    !*****************************************************************************************
    ! July, 2023
    !  - added the argument part
    !  - added allocation for vt(:), when part = 2
    implicit none
    type(conditional_ts), intent(out) :: cts
    integer,  intent(in) :: n, nreg, part
    real(dp), intent(in) :: xreg(n, max(1, nreg))
    real(dp), intent(in) :: xstart(max(1, nreg))

    !--------------------------------------------------------------
    ! the regressors and the corresponding starting values
    !--------------------------------------------------------------
    cts%nreg = nreg
    call safe_allocate(cts%xstart, max(1, nreg))
    cts%xstart = xstart
    if (nreg > 0) then
       call safe_allocate(cts%xreg, n, nreg)
       cts%xreg = xreg
    end if

    !--------------------------------------------------------------
    !   w - the conditional parameter
    !  gw - the transformed parameter (for now, only used in part 2)
    !--------------------------------------------------------------
    call safe_allocate(cts%w, n)
    cts%w  = 0.d0
    if (part == 2) then
       call safe_allocate(cts%gw, n)
       cts%gw  = 0.d0
    end if

    !--------------------------------------------------------------
    ! Transformed time series
    ! gi1 = g12(y) - might be used to calculate the error term
    ! gi2 = gi2    - the AR term
    !--------------------------------------------------------------
    if (part == 1) then
       call safe_allocate(cts%gi1, n)
       cts%gi1  = 0.d0
    end if
    call safe_allocate(cts%gi2, n)
    cts%gi2  = 0.d0

    !--------------------------------------------------------------
    ! Linear predictor
    ! eta = gi1(gw)
    !--------------------------------------------------------------
    call safe_allocate(cts%eta, n)
    cts%eta = 0.d0

    !--------------------------------------------------------------
    ! The error term
    ! et = g13(yt, wt),      part 1
    !      g23(g13(yt, wt)), part 2
    !--------------------------------------------------------------
    call safe_allocate(cts%et, n)
    cts%et  = 0.d0
    return
  end subroutine allocate_conditional_ts

  subroutine allocate_model_part(model, order, xreg, xstart, xregar, nfix, alpha, flagsb, &
       fvbeta, flagsar, fvar, flagsma, fvma, d, part)
    !*****************************************************************************************
    ! Allocates each part of the model:
    !  - parameters : allocates the parameter vectors and sets the fixed values
    !  - time series: allocates the conditional time series and sets initial values
    ! The value of "n" must be set in the object "model" before calling this subroutine
    !*****************************************************************************************
    ! Input
    !  order  : a size 4 vector with (nreg, p, q, inf)
    !  xreg   : matrix of regressors
    !  xstart : starting values for the regressors (to be used when t < 1)
    !  xregar : 0 = xreg is included only in the intercept
    !           1 = xreg is also included in the AR part.
    !  fita   : 1 = fit alpha, 0 = alpha is known
    !  alpha  : if fixa = 1, the value of the parameter alpha.
    !           if fixa = 0, a dummy argument.
    !  fitb   : number of unknown parameters in the beta vector (fitb <= nreg)
    !  flagsb : the lags that must be fixed in beta
    !  fvbeta : the value of the fixed parameters in beta
    !  flagsar: the lags that must be fixed in phi
    !  fvar   : the value of the fixed parameters in phi
    !  flagsma: the lags that must be fixed in theta
    !  fvma   : the value of the fixed parameters in theta
    !  d      : if fixd = 1, the value of the parameter alpha.
    !           if fixd = 0, a dummy argument.
    !  part   : indicates which part of the model is being allocated (1 or 2)
    !
    ! Input / Output
    !   model: input, the model's configurations partially defined.
    !          output, updated model's configurations
    !*****************************************************************************************
    ! July, 2023
    !   - fixnu no longer exist in the model object so the code to set
    !     this variable was removed
    !   - removed the argument m
    !   - changed d == 0.0d0 to abs(d) < epsilon(1.d0)
    implicit none
    type(argsModel), intent(inout) :: model
    integer,  target, intent(in) :: order(4), nfix(5)
    integer,  pointer    :: nreg, p, q, inf
    integer,  pointer    :: fixa, fixb, fixar, fixma, fixd
    integer,  intent(in) :: part, xregar
    integer,  intent(in) :: flagsb(max(1, nfix(2)))
    integer,  intent(in) :: flagsar(max(1, nfix(3)))
    integer,  intent(in) :: flagsma(max(1, nfix(4)))
    real(dp), intent(in) :: xstart(max(1, order(1)))
    real(dp), intent(in) :: alpha, d
    real(dp), intent(in) :: fvbeta(max(1, nfix(2)))
    real(dp), intent(in) :: fvar(max(1, nfix(3)))
    real(dp), intent(in) :: fvma(max(1, nfix(4)))
    real(dp), intent(in) :: xreg(model%n, max(1, order(1)))

    nreg => order(1)
    p    => order(2)
    q    => order(3)
    inf  => order(4)
    fixa  => nfix(1)
    fixb  => nfix(2)
    fixar => nfix(3)
    fixma => nfix(4)
    fixd  => nfix(5)

    !----------------------------------------------
    ! alpha, beta, ar, ma and d
    !----------------------------------------------
    call allocate_parvec(model%pt(part)%alpha, 1, fixa, [1], [alpha])
    call allocate_parvec(model%pt(part)%beta, nreg, fixb, flagsb, fvbeta)
    call allocate_parvec(model%pt(part)%ar, p, fixar, flagsar, fvar)
    call allocate_parvec(model%pt(part)%ma, q, fixma, flagsma, fvma)
    call allocate_parvec(model%pt(part)%d, 1, fixd, [1], [d])

    !---------------------------------------------
    !  If d = 0 uses inf = q
    !---------------------------------------------
    model%pt(part)%inf = max(inf, q)
    if (abs(d) < epsilon(1.d0) .and. fixd == 1) model%pt(part)%inf = q

    !----------------------------------------------
    ! mu, eta and xreg
    !----------------------------------------------
    call allocate_conditional_ts(model%cts(part), model%n, nreg, xreg, xstart, part)
    model%cts(part)%xregar = xregar
    return
  end subroutine allocate_model_part

  subroutine allocate_model(model, n, order, y, gy, xreg1, xreg2, xstart, xregar, nfix, alpha, &
       flagsb, fvbeta, flagsar, fvar, flagsma, fvma, d)
    !*****************************************************************************************
    ! Allocates the complete model. For each part
    !  - parameters : allocates the parameter vectors and sets the fixed values
    !  - time series: allocates the conditional time series and sets initial values
    !
    ! Vectors and matrices related to the score vector and information matrix are allocated on
    ! demand by subroutines.
    !*****************************************************************************************
    ! Input
    !  n      : sample size
    !  order  : a 2 by 4 matrix with (nreg, p, q, inf)
    !  y      : the observed time series
    !  gy     : a matrix n by 2 with the transformed time series g11(y) and g12(y)
    !  xreg1  : a matrix of regressors corresponding to part 1
    !  xreg2  : a matrix of regressors corresponding to part 2
    !  xstart : matrix of maxval(nreg) rows by 2 columns with starting values for the regressors
    !          (to be used when t < 1)
    !  xregar : a vector of size 2.
    !           0 = xreg is included only in the intercept
    !           1 = xreg is also included in the AR part.
    !  nfix   : a 2 by 5 matrix with the number of fixed parameters
    !  alpha  : a vector of size 2.
    !           if fita(i) = 0, the value of the parameter alpha corresponding to part i.
    !           if fita(i) = 1, a dummy argument.
    !  flagsb : matrix with two columns with the lags that must be fixed in each beta
    !  fvbeta : matrix with two columns with the value of the fixed parameters in each beta
    !  flagsar: matrix with two columns with the lags that must be fixed in each phi
    !  fvar   : matrix with two columns with the value of the fixed parameters in each phi
    !  flagsma: matrix with two columns with the lags that must be fixed in each theta
    !  fvma   : matrix with two columns with the value of the fixed parameters in each theta
    !  d      : if fixd = 1, the value of the parameter alpha.
    !           if fixd = 0, a dummy argument.
    !
    ! Input / Output
    !  model: input, the model's configurations pontentialy partially defined.
    !         output, updated model's configurations
    !*****************************************************************************************
    ! July, 2023
    !  - removed the subroutine allocate_model1 and renamed allocate_model2 to
    !    allocate_model
    !  - merged some arguments to simplify the code
    !     example: the vectors flagar1 and flagar2 are now one matrix named flagar)
    !  - changed the code accordingly
    !
    ! Last revision: March, 2025
    !  - now gy has 2 columns
    !  - replaced nreg, p, q and inf with the matrix order
    !  - replaced fita, fitb, fitar, fitma, fitd with nfix
    !    changed all subroutines accordingly
    !  - added pointers
    implicit none
    type(argsModel), intent(inout) :: model
    integer,  intent(in) :: order(2,4), nfix(2,5)
    integer,  intent(in) :: n, xregar(2)
    integer,  intent(in) :: flagsb(2, maxval([1, nfix(:, 2)]))
    integer,  intent(in) :: flagsar(2, maxval([1, nfix(:, 3)]))
    integer,  intent(in) :: flagsma(2, maxval([1, nfix(:, 4)]))
    real(dp), intent(in) :: alpha(2), d(2)
    real(dp), intent(in) :: fvbeta(2, maxval([1, nfix(:, 2)]))
    real(dp), intent(in) :: fvar(2, maxval([1, nfix(:, 3)]))
    real(dp), intent(in) :: fvma(2, maxval([1, nfix(:, 4)]))
    real(dp), intent(in) :: xreg1(n, max(1, order(1,1)))
    real(dp), intent(in) :: xreg2(n, max(1, order(2,1)))
    real(dp), intent(in) :: xstart(2, maxval([1, order(:,1)]))
    real(dp), intent(in) :: y(n), gy(n,2)
    real(dp), allocatable :: xreg(:,:)
    integer:: ib, ip, it, part

    !--------------------------------------------------------------
    ! Allocating y
    !--------------------------------------------------------------
    model%n = n
    call safe_allocate(model%y, n)
    model%y = y

    do part = 1,2
       ib = max(1, nfix(part,2))
       ip = max(1, nfix(part,3))
       it = max(1, nfix(part,4))

       !--------------------------------------------------------------
       ! Allocating each part
       !--------------------------------------------------------------
       call safe_allocate(xreg, n, max(1, order(part,1)))
       if(part == 1) then
          xreg = xreg1
       else
          xreg = xreg2
       end if
       call allocate_model_part(model, order(part,:), xreg, xstart(part,1:max(1, order(part,1))), &
            xregar(part), nfix(part,:), alpha(part), flagsb(part, 1:ib), fvbeta(part, 1:ib), &
            flagsar(part, 1:ip), fvar(part, 1:ip), flagsma(part, 1:it), fvma(part, 1:it), d(part), part)

       !--------------------------------------------------------------
       ! setting the number of parameters for each part of the model
       !--------------------------------------------------------------
       model%pt(part)%npar = 2 + sum(order(part, 1:3)) - sum(nfix(part,:))
    end do

    !--------------------------------------------------------------
    ! update g11(y)) and g12(y)
    !--------------------------------------------------------------
    model%cts(1)%gi1 = gy(:,1)
    model%cts(1)%gi2 = gy(:,2)

    !--------------------------------------------------------------
    ! update skip
    !--------------------------------------------------------------
    model%pt(1:2)%skip = 0
    if (model%pt(1)%npar == 0) model%pt(1)%skip = 1
    if (model%pt(2)%npar == 0) model%pt(2)%skip = 1
    return
  end subroutine allocate_model

  subroutine allocate_deta(deta, fita, fitb, fitar, fitma, fitd, n)
    !*****************************************************************************************
    ! Allocates deta/dalpha, deta/dbeta, deta/dphi, deta/dtheta and deta/dd
    !*****************************************************************************************
    ! Input
    !   fita : 1 = fit alpha, 0 = alpha is known
    !   fitb : number of unknown parameters in the beta vector
    !   fitar: number of unknown parameters in the phi (AR) vector
    !   fitma: number of unknown parameters in the theta (MA) vector
    !   fitd : 1 = fit d, 0 = d is known
    !   n    : sample size
    !
    ! Output
    !   deta: allocated vectors/matrices to save the values of deta/dgamma
    !*****************************************************************************************
    ! Last revision: July, 2023
    implicit none
    type(deta_d), intent(out) :: deta
    integer, intent(in) :: n, fita, fitb, fitar, fitma, fitd

    if (fita == 1) call safe_allocate(deta%dalpha, n, 1)
    if (fitb > 0)  call safe_allocate(deta%dbeta, n, fitb)
    if (fitar > 0) call safe_allocate(deta%dphi, n, fitar)
    if (fitma > 0) call safe_allocate(deta%dtheta, n, fitma)
    if (fitd == 1) call safe_allocate(deta%dd, n, 1)
    return
  end subroutine allocate_deta


  subroutine allocate_Us(U, fita, fitb, fitar, fitma, fitd)
    !*****************************************************************************************
    ! Allocates U(alpha), U(beta), U(phi), U(theta) and U(d)
    !*****************************************************************************************
    ! Input
    !   fita : 1 = fit alpha, 0 = alpha is known
    !   fitb : number of unknown parameters in the beta vector
    !   fitar: number of unknown parameters in the phi (AR) vector
    !   fitma: number of unknown parameters in the theta (MA) vector
    !   fitd : 1 = fit d, 0 = d is known
    !
    ! Output
    !   U: allocated vectors to save the values of the score vector
    !*****************************************************************************************
    ! Last revision: July, 2023
    !  - removed the argument fitnu and changed the code accordingly
    implicit none
    type(score), intent(inout) :: U
    integer, intent(in) :: fita, fitb, fitar, fitma, fitd

    if (fita == 1) call safe_allocate(U%Ualpha, 1)
    if (fitb > 0)  call safe_allocate(U%Ubeta, fitb)
    if (fitar > 0) call safe_allocate(U%Uphi, fitar)
    if (fitma > 0) call safe_allocate(U%Utheta, fitma)
    if (fitd == 1) call safe_allocate(U%Ud, 1)
    return
  end subroutine allocate_Us


  subroutine allocate_SI(model, SI)
    !*****************************************************************************************
    ! Allocates vectors and matrices related to the score vector and information matrix:
    !   T1, T2, h1, h2, deta1/drho, deta2/drho, deta1/dlambda and U
    !
    ! Input
    !   model: the model's configurations
    !
    ! Output
    !  SI: allocated vectors and matrices realted to the score vector and information matrix
    !*****************************************************************************************
    ! Last revision: July, 2023
    implicit none
    type(argsModel), intent(in) :: model
    type(argsSI), intent(out) :: SI
    integer :: n1, n2

    n1 = min(model%n, model%n * (1 - model%pt(1)%skip) + 1)
    n2 = min(model%n, model%n * (1 - model%pt(2)%skip) + 1)

    !--------------------------------------------------------------------
    ! T and h - for compatibility these vectors always must be allocated
    !--------------------------------------------------------------------
    call safe_allocate(SI%T(1)%z, n1)
    call safe_allocate(SI%h(1)%z, n1)
    call safe_allocate(SI%T(2)%z, n2)
    call safe_allocate(SI%h(2)%z, n2)

    if (model%pt(1)%npar > 0) then
       !--------------------------------------------------------------------
       ! deta1_drho
       !--------------------------------------------------------------------
       call allocate_deta(SI%deta(1, 1), model%pt(1)%alpha%fit, model%pt(1)%beta%fit, &
            model%pt(1)%ar%fit, model%pt(1)%ma%fit, model%pt(1)%d%fit, model%n)

       !--------------------------------------------------------------------
       ! Urho
       !--------------------------------------------------------------------
       if (model%sco == 1) call allocate_Us(SI%U(1), model%pt(1)%alpha%fit, &
            model%pt(1)%beta%fit, model%pt(1)%ar%fit, model%pt(1)%ma%fit, model%pt(1)%d%fit)
    end if

    if (model%pt(2)%npar == 0) return

    !--------------------------------------------------------------------
    ! deta2_drho
    !--------------------------------------------------------------------
    if (model%pt(1)%npar > 0) call allocate_deta(SI%deta(2, 1), model%pt(1)%alpha%fit, &
         model%pt(1)%beta%fit, model%pt(1)%ar%fit, model%pt(1)%ma%fit, model%pt(1)%d%fit, model%n)

    !--------------------------------------------------------------------
    ! deta2_dlambda
    !--------------------------------------------------------------------
    call allocate_deta(SI%deta(2, 2), model%pt(2)%alpha%fit, &
         model%pt(2)%beta%fit, model%pt(2)%ar%fit, model%pt(2)%ma%fit, model%pt(2)%d%fit, model%n)

    !--------------------------------------------------------------------
    ! Ulambda
    !--------------------------------------------------------------------
    if (model%sco == 1) call allocate_Us(SI%U(2), model%pt(2)%alpha%fit, &
         model%pt(2)%beta%fit, model%pt(2)%ar%fit, model%pt(2)%ma%fit, model%pt(2)%d%fit)
    return
  end subroutine allocate_SI

end module alloc
