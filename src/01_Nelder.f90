MODULE Nelder_mead
  !-----------------------------------------------------------------------------------
  ! Modifyed by: Taiane Schaedler Prass
  ! Porto Alegre, 8/2019
  !
  ! import the double precision kind from the main module that contains the 
  ! definition of the special type variables
  !-----------------------------------------------------------------------------------
  use main_mod, only : argsModel

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: minim

  integer, parameter :: dp = kind(1.d0)  

CONTAINS

!!$  function wname(act, i1, db1, db2) result(output)   
!!$    implicit none
!!$    character(len = *), intent(in) :: act
!!$    integer, intent(in) :: i1
!!$    real(dp), intent(in) :: db1, db2
!!$    character(len = 15) :: out1, out2, out3
!!$    character(len = 100) :: output
!!$
!!$    write(out1,fmt = "(i7)") i1
!!$    write(out2, fmt = "(f12.4)") db1
!!$    write(out3,fmt = "(f12.4)") db2
!!$
!!$    output = trim(adjustl(act))//" / "//trim(adjustl(out1))//" / "//&
!!$         trim(adjustl(out2))//" / "//trim(adjustl(out3))
!!$    return
!!$  end function wname


  SUBROUTINE minim(p, step, nop, func, maxfn, iprint, stopcr, functn, &
       ifault, neval, model)
    !
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
    !     Modified by taiane schaedler prass  (08/2019).
    !
    !     This subroutine was modified to have as imput the new type of variable "model"
    !     Which is used to pass extra arguments for the log-likelihood function:
    !
    !     Functn  = input, name of the user's subroutine - arguments (model,nop,p,func)
    !               Which returns the function value for a given set of
    !               Parameter values in array p.
    !
    !     Neval = output, number of function evaluations performed by the algorithm
    !   
    !     Model = input, user defined variable with extra information to be passed to functn 
    !**********************************************************************************************      

    type(argsModel)            :: model 
    INTEGER, INTENT(IN)        :: nop, maxfn, iprint
    INTEGER, INTENT(OUT)       :: ifault, neval
    REAL (dp), INTENT(IN)      :: stopcr
    REAL (dp), INTENT(IN OUT)  :: p(nop), step(nop)
    REAL (dp), INTENT(OUT)     :: func

    INTERFACE
       SUBROUTINE functn(model, nop, p, func)
         import :: dp, argsModel 
         IMPLICIT NONE
         type(argsModel), intent(inout) :: model
         integer, intent(in) :: nop
         REAL (dp), INTENT(IN)  :: p(nop)
         REAL (dp), INTENT(OUT) :: func
       END SUBROUTINE functn
    END INTERFACE

    !     Local variables

    REAL (dp)   :: g(nop+1,nop), h(nop+1), pbar(nop), pstar(nop), pstst(nop), &
         hstst, hmin, hmean, hstd, hstar, hmax, savemn

    REAL (dp), PARAMETER :: zero = 0.d0, half = 0.5d0, one = 1.d0, two = 2.d0
    INTEGER     :: i, imax, imin, irow, j, nap,  np1, iflag

    !     A = REFLECTION COEFFICIENT, B = CONTRACTION COEFFICIENT, AND
    !     C = EXPANSION COEFFICIENT.

    REAL (dp), PARAMETER :: a = 1.d0, b = 0.5d0, c = 2.d0, eps = 0.001d0

    character(len = 100) :: output

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
    IF (nop < 1) ifault = 2
    IF (ifault /= 0) RETURN

    !     SET NAP = NO. OF PARAMETERS TO BE VARIED, I.E. WITH STEP /= 0

    nap = count(step /= zero)
    neval = 0
    np1 = nap + 1
    iflag = 0

    !     IF NAP = 0 EVALUATE FUNCTION AT THE STARTING POINT AND RETURN

    IF (nap <= 0 .or. maxfn <= 1)then
       CALL functn(model,nop,p,func)
       RETURN
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
    DO i = 1, nop
       IF (step(i) /= zero) THEN
          g(irow,:) = p
          g(irow,i) = p(i) + step(i)
          irow = irow + 1
       END IF
    END DO

    DO i = 1, np1
       p = g(i,:)
       CALL functn(model,nop,p,h(i))
       neval = neval + 1
    END DO

    hmin = h(1)
    if(iprint > 0) then
       DO i = 1, np1
          hmin = min(hmin, h(i))
          !output = wname("  BUILD", i, h(i), hmin)
          !call labelpr(output, -1)
          call dblepr("  BUILD", -1, (/dble(i), h(i), hmin/), 3)
       end DO
    end if

    !     START OF MAIN CYCLE.
    !     Do while hstd > stopcr .AND. neval <= maxfn
    Main_loop: DO      

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
       DO i = 2, np1
          IF (h(i) > hmax) THEN
             imax = i
             hmax = h(i)
          ELSE
             IF (h(i) < hmin) THEN
                imin = i
                hmin = h(i)
             END IF
          END IF
       END DO

       !------------------------------------------
       !  Check to see if minimum reached.
       !------------------------------------------
       !     CALCULATE MEAN & STANDARD DEVIATION OF FUNCTION VALUES FOR THE
       !     CURRENT SIMPLEX.
       hmean = SUM(h)/np1
       hstd = SUM((h - hmean)** 2)
       hstd = SQRT(hstd/np1)

       !     IF THE RMS > STOPCR, SET IFLAG & LOOP TO ZERO AND GO TO THE
       !     START OF THE MAIN CYCLE AGAIN.
       if(neval > maxfn) exit main_loop

       IF (hstd > stopcr) iflag = 0
       !     CONVERGENCE CRITERION SATISFIED.
       !     IF IFLAG = 0, SET IFLAG & SAVE HMEAN.
       !     IF IFLAG = 1 & CHANGE IN HMEAN <= STOPCR THEN SEARCH IS COMPLETE.
       IF (iflag == 0 .OR. ABS(savemn-hmean) >= stopcr) THEN
          iflag = 1
          savemn = hmean
       ELSE
          EXIT Main_loop
       END IF

       !     FIND THE CENTROID OF THE VERTICES OTHER THAN P(IMAX)

       pbar = zero
       DO i = 1, np1
          IF (i /= imax) THEN
             pbar = pbar + g(i,:)
          END IF
       END DO
       pbar = pbar / nap

       !------------------------------------------------
       !  Reflection through the centroid.
       !------------------------------------------------
       !     REFLECT MAXIMUM THROUGH PBAR TO PSTAR,
       !     HSTAR = FUNCTION VALUE AT PSTAR.

       pstar = pbar +  a * (pbar -  g(imax,:))
       CALL functn(model,nop,pstar,hstar)
       neval = neval + 1

       IF (hstar < hmin) THEN         

          !-------------------------------------------------
          !  Successful reflection, so extension.
          !-------------------------------------------------
          !     IF HSTAR < HMIN, REFLECT PBAR THROUGH PSTAR,
          !     HSTST = FUNCTION VALUE AT PSTST.
          pstst = pbar + c * (pstar - pbar) 
          CALL functn(model,nop,pstst,hstst)
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
             IF (iprint > 0) THEN
                IF (MOD(neval,iprint) == 0) then
                   !output = wname("  EXTENSION", neval, hstst, hmin)
                   !call labelpr(output, -1)
                   call dblepr("  EXTENSION", -1, (/dble(neval), hstst, hmin/), 3)
                end if
             END IF
          else   
             ! REPLACE MAXIMUM POINT BY PSTAR & H(IMAX) BY HSTAR.
             g(imax,:) = pstar
             h(imax) = hstar   
          END IF
          goto 250
       end IF

       !     HSTAR IS NOT < HMIN 
       !     TEST WHETHER IT IS < FUNCTION VALUE AT SOME POINT OTHER THAN
       !     P(IMAX).   IF IT IS REPLACE P(IMAX) BY PSTAR & HMAX BY HSTAR.
       DO i = 1, np1
          IF (i /= imax) THEN
             IF (hstar < h(i)) THEN  
                !-------------------------------------------------
                !  pstar is better than the second best point.
                !  Accept reflection and restart the loop.
                !-------------------------------------------------    
                g(imax,:) = pstar
                h(imax) = hstar
                IF (iprint > 0) THEN
                   IF (MOD(neval,iprint) == 0) then
                      !output = wname("  REFLECTION", neval, hstar, hmin)
                      !call labelpr(output, -1)
                      call dblepr("  REFLECTION", -1, (/dble(neval), hstar, hmin/), 3)
                   end if
                END IF
                GO TO 250
             END IF
          END IF
       END DO

       !     HSTAR > ALL FUNCTION VALUES EXCEPT POSSIBLY HMAX.
       !     IF HSTAR <= HMAX, REPLACE P(IMAX) BY PSTAR & HMAX BY HSTAR.
       !     OTHERWHISE, THE OLD VALUE OF P(IMAX) IS USED TO
       !     OBTAIN THE REFLECTION

       j = 0 ! (Inside contraction)
       IF (hstar < hmax) THEN
          g(imax,:) = pstar
          hmax = hstar
          h(imax) = hstar
          j = 1   ! (Outside contraction)
       END IF

       !-----------------------------------------------------
       !  Contraction using the best of the two points
       !  Ouside-Contraction
       !-----------------------------------------------------
       !     CONTRACTED STEP TO THE POINT PSTST,
       !     HSTST = FUNCTION VALUE AT PSTST.
       pstst =  pbar +  b * (g(imax,:) - pbar)
       CALL functn(model,nop,pstst,hstst)
       neval = neval + 1

       if(hstst <= hmax) then
          !-----------------------------
          !  Retain contraction.
          !-----------------------------
          !     IF HSTST < HMAX REPLACE P(IMAX) BY PSTST & HMAX BY HSTST.
          g(imax,:) = pstst
          h(imax) = hstst               

          IF (iprint > 0) THEN
             IF (MOD(neval,iprint) == 0) then                
                !-----------------------------------------------------
                !  Ouside-Contraction
                !-----------------------------------------------------
                !if(j == 1) output = wname("  OUT-REDUCTION", neval, hstst, hmin)
                if(j == 1) call dblepr("  OUT-REDUCTION", -1, (/dble(neval), hstst, hmin/), 3)
                !-----------------------------------------------------
                !  Inside-contraction
                !-----------------------------------------------------
                if(j == 0) call dblepr("  In-REDUCTION", -1, (/dble(neval), hstst, hmin/), 3)
                !if(j == 0)output = wname("  In-REDUCTION", neval, hstst, hmin)
                !call labelpr(output, -1)
             end if
          END IF
          goto 250
       end if

       !     HSTST > HMAX.
       !--------------------------------------
       !  Contract the whole simplex.
       !--------------------------------------
       !     SHRINK THE SIMPLEX BY REPLACING EACH POINT, OTHER THAN THE CURRENT
       !     MINIMUM, BY A POINT MID-WAY BETWEEN ITS CURRENT POSITION AND THE
       !     MINIMUM.

       DO i = 1, np1
          IF (i /= imin) THEN
             DO j = 1, nop
                IF (step(j) /= zero) g(i,j) = (g(i,j) + g(imin,j))*half
                p(j) = g(i,j)
             END DO
             CALL functn(model,nop,p,h(i))
             neval = neval + 1
             IF (iprint > 0) THEN
                IF (mod(neval,iprint) == 0) then
                   !output = wname("  SHRINK", neval, h(i), hmin)
                   !call labelpr(output, -1)
                   call dblepr("  SHRINK", -1, (/dble(neval), h(i), hmin/), 3)
                end if
             END IF
          END IF
       END DO

    END DO Main_loop

    !--------------------------------------------------------------
    !  Minimum of last simplex
    !--------------------------------------------------------------
    p = g(imin,:)
    func = h(imin)


    !     TEST WHETHER THE NO. OF FUNCTION VALUES ALLOWED, maxfn, HAS BEEN
    !     OVERRUN; IF SO, EXIT WITH IFAULT = 1.

    IF (neval >= maxfn) THEN
       ifault = 1     
       IF (iprint < 0) RETURN
       call intpr1(" No. of function evaluations > ", -1,  maxfn)
       call dblepr1(" RMS of function values of last simplex =", -1, hstd)
       call dblepr(" Minimum in the last simplex", -1, p, nop)
       call dblepr1(" Function value at minimum =", -1, func)            
       RETURN
    END IF

    !     CONVERGENCE CRITERION SATISFIED.
    IF (iprint >= 0) THEN
       call labelpr(" EVIDENCE OF CONVERGENCE", -1)
       call intpr1(" No. of function evaluations = ", -1,  neval)
       call dblepr1(" RMS of function values of last simplex =", -1, hstd)
       call dblepr(" Minimum at", -1, p, nop)
       call dblepr1(" Function value at minimum =", -1, func)    
    END IF

    RETURN
  END SUBROUTINE minim

END MODULE Nelder_mead
