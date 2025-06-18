module Nelder_mead
  use main_mod, only : optimFunc, argsModel, dp
  !-------------------------------------------------------------------------------------------
  ! Modifyed by: Taiane Schaedler Prass
  ! Porto Alegre, 8/2019
  !
  ! depends on the main module that contains the definition of the special type variables
  !-------------------------------------------------------------------------------------------
  ! Last revision: February 2024
  !  - no changes from previous version
  !-------------------------------------------------------------------------------------------
  implicit none
  private
  public :: minim

contains

  subroutine minim(p, step, nop, func, maxfn, iprint, stopcr, loglik, ifault, neval, model)
    !********************************************************************************************
    !     A program for function minimization using the simplex method.
    !
    !  Reference:
    !
    !    John Nelder, Roger Mead,
    !    A simplex method for function minimization,
    !    Computer Journal,
    !    Volume 7, 1965, pages 308-313.
    !
    !    R ONeill,
    !    Algorithm AS 47:
    !    Function Minimization Using a Simplex Procedure,
    !    Applied Statistics,
    !    Volume 20, Number 3, 1971, pages 338-345.
    !
    !---------------
    !     Programmed by D.E. Shaw,
    !     Csiro, Division of Mathematics & Statistics
    !     P.o. Box 218, Lindfield, N.S.W. 2070
    !---------------
    !     With amendments by R.W.M. Wedderburn
    !     Rothamsted Experimental Station
    !     Harpenden, Hertfordshire, England
    !---------------
    !     Further amended by Alan Miller
    !     Csiro Division of Mathematical & Information Sciences
    !     Private bag 10, Clayton, Vic. 3169
    !---------------
    !     Fortran 90 conversion by Alan Miller, june 1995
    !     Alan.miller @ vic.cmis.csiro.au
    !     Latest revision - 5 december 1999
    !---------------
    !     Further modified by
    !     Taiane Schaedler Prass, april 2020
    !     Statistics Department, PPGEst, UFRGS
    !     taiane.prass@ufrgs.br
    !     (for detail see http://www.scholarpedia.org/article/Nelder-Mead_algorithm)
    !-----------------
    !     Arguments:-
    !     P()     = input, starting values of parameters
    !               Output, final values of parameters
    !     Step()  = input, initial step sizes
    !     Nop     = input, no. Of parameters, incl. Any to be held fixed
    !     Func    = output, the function value corresponding to the final
    !                 Parameter values.
    !     Maxfn     = input, the maximum no. Of function evaluations allowed.
    !               Say, 20 times the number of parameters, nop.
    !     Iprint  = input, print control parameter
    !                 < 0 No printing
    !                 = 0 Printing of parameter values and the function
    !                     Value after initial evidence of convergence.
    !                 > 0 As for iprint = 0 plus progress reports after every
    !                     Iprint evaluations, plus printing for the initial simplex.
    !     Stopcr  = input, stopping criterion.
    !               The criterion is applied to the standard deviation of
    !               The values of func at the points of the simplex.
    !     Iquad   = input, = 1 if fitting of a quadratic surface is required
    !                      = 0 If not
    !               N.b. The fitting of a quadratic surface is strongly
    !               Recommended, provided that the fitted function is
    !               Continuous in the vicinity of the minimum.   It is often
    !               A good indicator of whether a premature termination of
    !               The search has occurred.
    !     Simp    = input, criterion for EXTENSION the simplex to overcome
    !               Rounding errors before fitting the quadratic surface.
    !               The simplex is expanded so that the function values at
    !               The points of the simplex exceed those at the supposed
    !               Minimum by at least an amount simp.
    !     Var()   = output, contains the diagonal elements of the inverse of
    !               The information matrix.
    !     Functn  = input, name of the user's subroutine - arguments (p,func)
    !               Which returns the function value for a given set of
    !               Parameter values in array p.
    !****     Functn must be declared external in the calling program.
    !     Ifault  = output, = 0 for successful termination
    !                 = 1 If maximum no. Of function evaluations exceeded
    !                 = 2 If nop < 1 or stopcr <= 0
    !
    !     N.b. P, step and var (if iquad = 1) must have dimension at least nop
    !          In the calling program.
    !
    !********************************************************************************************
    !
    !  Modified by taiane schaedler prass  (08/2019).
    !
    !  This subroutine was modified to have as imput the new type of variable "model"
    !  which is used to pass extra arguments for the log-likelihood function:
    !
    !  Functn  = input, name of the user's subroutine - arguments (model,nop,p,func)
    !            Which returns the function value for a given set of
    !            Parameter values in array p.
    !  Neval = output, number of function evaluations performed by the algorithm
    !  Model = input, user defined variable with extra information to be passed to functn
    !
    !  Last update: March, 2025
    !    - added the optimFunc user defined type so the interface can be removed
    !**********************************************************************************************

    type(optimFunc), intent(inout)    :: loglik
    type(argsmodel), intent(inout) :: model
    integer,   intent(in)          :: nop, maxfn, iprint
    integer,   intent(out)         :: ifault, neval
    real (dp), intent(in)          :: stopcr
    real (dp), intent(inout)       :: p(nop), step(nop)
    real (dp), intent(out)         :: func

    !     local variables

    real (dp)   :: g(nop+1,nop), h(nop+1), pbar(nop), pstar(nop), pstst(nop), &
         hstst, hmin, hmean, hstd, hstar, hmax, savemn

    real (dp), parameter :: zero = 0.d0, half = 0.5d0, one = 1.d0, two = 2.d0
    integer     :: i, imax, imin, irow, j, nap,  np1, iflag

    !     A = REFLECTION COEFFICIENT, B = CONTRACTION COEFFICIENT, AND
    !     C = EXPANSION COEFFICIENT.

    real (dp), parameter :: a = 1.d0, b = 0.5d0, c = 2.d0, eps = 0.001d0

    savemn = 0.d0 ! initialization to avoid compiler warning

    !------------------------------------------------------------------
    !     IF PROGRESS REPORTS HAVE BEEN REQUESTED, PRINT HEADING
    !------------------------------------------------------------------

    if(iprint >= 0) then
       call labelpr("------------------------------------------------", -1)
       call labelpr("  Nelder-Mead direct search function minimizer  ", -1)
       call labelpr("------------------------------------------------", -1)
    end if
    if(iprint > 0) then
       call labelpr("Progress Report every 'iprint' function evaluations", -1)
       call intpr1("iprint = ", -1, iprint)
       call labelpr(" ACTION / NEVAL / FUNC.VALUE. / LOWER", -1)
    end if

    !     CHECK INPUT ARGUMENTS

    ifault = 0
    if (stopcr <= 0.d0) ifault = 2
    if (nop < 1) ifault = 2
    if (ifault /= 0) return

    !     SET NAP = NO. OF PARAMETERS TO BE VARIED, I.E. WITH STEP /= 0

    nap = count(step /= zero)
    neval = 0
    np1 = nap + 1
    iflag = 0

    !     IF NAP = 0 EVALUATE FUNCTION AT THE STARTING POINT AND RETURN

    if (nap <= 0 .or. maxfn <= 1)then
       call loglik%functn(model,nop,p,func)
       return
    end if

    !-----------------------------------------
    !  Initial or restarted loop.
    !-----------------------------------------

    !     SET UP THE INITIAL SIMPLEX

    g = 0.d0
    h = 0.d0

    !20  continue

    !
    !  Define the initial simplex.
    !
    g(1,:) = p
    irow = 2
    do i = 1, nop
       if (step(i) /= zero) then
          g(irow,:) = p
          g(irow,i) = p(i) + step(i)
          irow = irow + 1
       end if
    end do

    do i = 1, np1
       p = g(i,:)
       call loglik%functn(model,nop,p,h(i))
       neval = neval + 1
    end do

    hmin = h(1)
    if(iprint > 0) then
       do i = 1, np1
          hmin = min(hmin, h(i))
          call dblepr("  build", -1, [dble(i), h(i), hmin], 3)
       end do
    end if

    !     START OF MAIN CYCLE.
    !     Do while hstd > stopcr .AND. neval <= maxfn
    main_loop: do

250    continue

       !     FIND MAX. & MIN. VALUES FOR CURRENT SIMPLEX (HMAX & HMIN).
       !
       !  Find highest and lowest function values.
       !  hmax = h(imax) indicates the vertex of the
       !  simplex to be replaced.
       !
       imax = 1
       imin = 1
       hmax = h(1)
       hmin = h(1)
       do i = 2, np1
          if (h(i) > hmax) then
             imax = i
             hmax = h(i)
          else
             if (h(i) < hmin) then
                imin = i
                hmin = h(i)
             end if
          end if
       end do

       !------------------------------------------
       !  Check to see if minimum reached.
       !------------------------------------------
       !     CALCULATE MEAN & STANDARD DEVIATION OF FUNCTION VALUES FOR THE
       !     CURRENT SIMPLEX.
       hmean = sum(h)/np1
       hstd = sum((h - hmean)** 2)
       hstd = sqrt(hstd/np1)

       !     IF THE RMS > STOPCR, SET IFLAG & LOOP TO ZERO AND GO TO THE
       !     START OF THE MAIN CYCLE AGAIN.
       if(neval > maxfn) exit main_loop

       if (hstd > stopcr) iflag = 0
       !     CONVERGENCE CRITERION SATISFIED.
       !     IF IFLAG = 0, SET IFLAG & SAVE HMEAN.
       !     IF IFLAG = 1 & CHANGE IN HMEAN <= STOPCR THEN SEARCH IS COMPLETE.
       if (iflag == 0 .or. abs(savemn-hmean) >= stopcr) then
          iflag = 1
          savemn = hmean
       else
          exit main_loop
       end if

       !     FIND THE CENTROID OF THE VERTICES OTHER THAN P(IMAX)

       pbar = zero
       do i = 1, np1
          if (i /= imax) then
             pbar = pbar + g(i,:)
          end if
       end do
       pbar = pbar / nap

       !------------------------------------------------
       !  Reflection through the centroid.
       !------------------------------------------------
       !     REFLECT MAXIMUM THROUGH PBAR TO PSTAR,
       !     HSTAR = FUNCTION VALUE AT PSTAR.

       pstar = pbar +  a * (pbar -  g(imax,:))
       call loglik%functn(model,nop,pstar,hstar)
       neval = neval + 1

       if (hstar < hmin) then

          !-------------------------------------------------
          !  Successful reflection, so extension.
          !-------------------------------------------------
          !     IF HSTAR < HMIN, REFLECT PBAR THROUGH PSTAR,
          !     HSTST = FUNCTION VALUE AT PSTST.
          pstst = pbar + c * (pstar - pbar)
          call loglik%functn(model,nop,pstst,hstst)
          neval = neval + 1

          !----------------------------------------------
          !  Retain extension?
          !----------------------------------------------
          !  Taiane Schaedler Prass - UFRGS
          !  Porto Alegre, 08/04/2020
          !  Replaced the condiction HSTST < HMIN (original code)
          !                       by HSTST < hstar
          !-----------------------------------------------
          if(hstst < hstar) then
             !     IF HSTST < HSTAR REPLACE CURRENT MAXIMUM POINT BY PSTST AND
             !     HMAX BY HSTST, THEN TEST FOR CONVERGENCE.
             g(imax,:) = pstst
             h(imax) = hstst
             if (iprint > 0) then
                if (mod(neval,iprint) == 0) then
                   call dblepr("  EXTENSION", -1, [dble(neval), hstst, hmin], 3)
                end if
             end if
          else
             ! REPLACE MAXIMUM POINT BY PSTAR & H(IMAX) BY HSTAR.
             g(imax,:) = pstar
             h(imax) = hstar
          end if
          goto 250
       end if

       !     HSTAR IS NOT < HMIN
       !     TEST WHETHER IT IS < FUNCTION VALUE AT SOME POINT OTHER THAN
       !     P(IMAX).   IF IT IS REPLACE P(IMAX) BY PSTAR & HMAX BY HSTAR.
       do i = 1, np1
          if (i /= imax) then
             if (hstar < h(i)) then
                !-------------------------------------------------
                !  pstar is better than the second best point.
                !  Accept reflection and restart the loop.
                !-------------------------------------------------
                g(imax,:) = pstar
                h(imax) = hstar
                if (iprint > 0) then
                   if (mod(neval,iprint) == 0) then
                      call dblepr("  REFLECTION", -1, [dble(neval), hstar, hmin], 3)
                   end if
                end if
                go to 250
             end if
          end if
       end do

       !     HSTAR > ALL FUNCTION VALUES EXCEPT POSSIBLY HMAX.
       !     IF HSTAR <= HMAX, REPLACE P(IMAX) BY PSTAR & HMAX BY HSTAR.
       !     OTHERWHISE, THE OLD VALUE OF P(IMAX) IS USED TO
       !     OBTAIN THE REFLECTION

       j = 0 ! (Inside contraction)
       if (hstar < hmax) then
          g(imax,:) = pstar
          hmax = hstar
          h(imax) = hstar
          j = 1   ! (Outside contraction)
       end if

       !-----------------------------------------------------
       !  Contraction using the best of the two points
       !  Ouside-Contraction
       !-----------------------------------------------------
       !     CONTRACTED STEP TO THE POINT PSTST,
       !     HSTST = FUNCTION VALUE AT PSTST.
       pstst =  pbar +  b * (g(imax,:) - pbar)
       call loglik%functn(model,nop,pstst,hstst)
       neval = neval + 1

       if(hstst <= hmax) then
          !-----------------------------
          !  Retain contraction.
          !-----------------------------
          !     IF HSTST < HMAX REPLACE P(IMAX) BY PSTST & HMAX BY HSTST.
          g(imax,:) = pstst
          h(imax) = hstst

          if (iprint > 0) then
             if (mod(neval,iprint) == 0) then
                !-----------------------------------------------------
                !  Ouside-Contraction
                !-----------------------------------------------------
                if(j == 1) call dblepr("  OUT-REDUCTION", -1, [dble(neval), hstst, hmin], 3)
                !-----------------------------------------------------
                !  Inside-contraction
                !-----------------------------------------------------
                if(j == 0) call dblepr("  In-REDUCTION", -1, [dble(neval), hstst, hmin], 3)
             end if
          end if
          goto 250
       end if

       !     HSTST > HMAX.
       !--------------------------------------
       !  Contract the whole simplex.
       !--------------------------------------
       !     SHRINK THE SIMPLEX BY REPLACING EACH POINT, OTHER THAN THE CURRENT
       !     MINIMUM, BY A POINT MID-WAY BETWEEN ITS CURRENT POSITION AND THE
       !     MINIMUM.

       do i = 1, np1
          if (i /= imin) then
             do j = 1, nop
                if (step(j) /= zero) g(i,j) = (g(i,j) + g(imin,j))*half
                p(j) = g(i,j)
             end do
             call loglik%functn(model,nop,p,h(i))
             neval = neval + 1
             if (iprint > 0) then
                if (mod(neval,iprint) == 0) then
                   call dblepr("  SHRINK", -1, [dble(neval), h(i), hmin], 3)
                end if
             end if
          end if
       end do

    end do main_loop

    !--------------------------------------------------------------
    !  Minimum of last simplex
    !--------------------------------------------------------------
    p = g(imin,:)
    func = h(imin)

    !     TEST WHETHER THE NO. OF FUNCTION VALUES ALLOWED, maxfn, HAS BEEN
    !     OVERRUN; IF SO, EXIT WITH IFAULT = 1.

    if (neval >= maxfn) then
       ifault = 1
       if (iprint < 0) return
       call intpr1(" No. of function evaluations > ", -1,  maxfn)
       call dblepr1(" RMS of function values of last simplex =", -1, hstd)
       call dblepr(" Minimum in the last simplex", -1, p, nop)
       call dblepr1(" Function value at minimum =", -1, func)
       return
    end if

    !     CONVERGENCE CRITERION SATISFIED.
    if (iprint >= 0) then
       call labelpr(" EVIDENCE OF CONVERGENCE", -1)
       call intpr1(" No. of function evaluations = ", -1,  neval)
       call dblepr1(" RMS of function values of last simplex =", -1, hstd)
       call dblepr(" Minimum at", -1, p, nop)
       call dblepr1(" Function value at minimum =", -1, func)
    end if

    return
  end subroutine minim

end module nelder_mead
