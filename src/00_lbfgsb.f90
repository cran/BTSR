module lbfgsb
  use lpk
  !---------------------------------------------------------------------------
  !
  ! contains the suroutines of l-bfgs-b algorithm
  !
  ! original code: http://users.iems.northwestern.edu/~nocedal/lbfgsb.html
  !
  ! Converted to F90 by Taiane Schaedler Prass
  !-----------------------------------------------------------------------------
  implicit none
  private
  public :: setulb

  integer, parameter :: dp = kind(1.d0)

contains

  !c                                                                                      
  !c  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”        
  !c  or “3-clause license”)                                                              
  !c  Please read attached file License.txt                                               
  !c                                        
  !c===========   L-BFGS-B (version 3.0.  April 25, 2011  ===================
  !c
  !c     This is a modified version of L-BFGS-B. Minor changes in the updated 
  !c     code appear preceded by a line comment as follows 
  !c  
  !c     c-jlm-jn 
  !c
  !c     Major changes are described in the accompanying paper:
  !c
  !c         Jorge Nocedal and Jose Luis Morales, Remark on "Algorithm 778: 
  !c         L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained 
  !c         Optimization"  (2011). To appear in  ACM Transactions on 
  !c         Mathematical Software,
  !c
  !c     The paper describes an improvement and a correction to Algorithm 778. 
  !c     It is shown that the performance of the algorithm can be improved 
  !c     significantly by making a relatively simple modication to the subspace 
  !c     minimization phase. The correction concerns an error caused by the use 
  !c     of routine dpmeps to estimate machine precision. 
  !c
  !c     The total work space **wa** required by the new version is 
  !c 
  !c                  2*m*n + 11m*m + 5*n + 8*m 
  !c
  !c     the old version required 
  !c
  !c                  2*m*n + 12m*m + 4*n + 12*m 
  !c
  !c
  !c            J. Nocedal  Department of Electrical Engineering and
  !c                        Computer Science.
  !c                        Northwestern University. Evanston, IL. USA
  !c
  !c
  !c           J.L Morales  Departamento de Matematicas, 
  !c                        Instituto Tecnologico Autonomo de Mexico
  !c                        Mexico D.F. Mexico.
  !c
  !c                        March  2011    
  !c                                                 
  !c============================================================================= 
  subroutine setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, &
       task, iprint, csave, lsave, isave, dsave)
    implicit none
    character (len=60), intent(inout) :: task, csave
    logical, intent(inout) :: lsave(4)
    integer, intent(in) :: n,m, iprint
    integer, intent(inout) :: nbd(n), iwa(3*n), isave(44)
    real (dp), intent(in) ::  factr, pgtol
    real (dp), intent(inout) :: x(n), l(n), u(n), f, g(n)
    !c-jlm-jn		
    real (dp), intent(inout) :: wa(2*m*n + 5*n + 11*m*m + 8*m), dsave(29)       
    !
    !c      ************
    !c 
    !c      Subroutine setulb
    !c 
    !c      This subroutine partitions the working arrays wa and iwa, and 
    !c        then uses the limited memory BFGS method to solve the bound
    !c        constrained optimization problem by calling mainlb.
    !c        (The direct method will be used in the subspace minimization.)
    !c 
    !c      n is an integer variable.
    !c        On entry n is the dimension of the problem.
    !c        On exit n is unchanged.
    !c 
    !c      m is an integer variable.
    !c        On entry m is the maximum number of variable metric corrections
    !c          used to define the limited memory matrix.
    !c        On exit m is unchanged.
    !c 
    !c      x is a double precision array of dimension n.
    !c        On entry x is an approximation to the solution.
    !c        On exit x is the current approximation.
    !c 
    !c      l is a double precision array of dimension n.
    !c        On entry l is the lower bound on x.
    !c        On exit l is unchanged.
    !c 
    !c      u is a double precision array of dimension n.
    !c        On entry u is the upper bound on x.
    !c        On exit u is unchanged.
    !c 
    !c      nbd is an integer array of dimension n.
    !c        On entry nbd represents the type of bounds imposed on the
    !c          variables, and must be specified as follows:
    !c          nbd(i)=0 if x(i) is unbounded,
    !c                 1 if x(i) has only a lower bound,
    !c                 2 if x(i) has both lower and upper bounds, and
    !c                 3 if x(i) has only an upper bound.
    !c        On exit nbd is unchanged.
    !c 
    !c      f is a double precision variable.
    !c        On first entry f is unspecified.
    !c        On final exit f is the value of the function at x.
    !c 
    !c      g is a double precision array of dimension n.
    !c        On first entry g is unspecified.
    !c        On final exit g is the value of the gradient at x.
    !c 
    !c      factr is a double precision variable.
    !c        On entry factr >= 0 is specified by the user.  The iteration
    !c          will stop when
    !c 
    !c          (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
    !c 
    !c          where epsmch is the machine precision, which is automatically
    !c          generated by the code. Typical values for factr: 1.d+12 for
    !c          low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
    !c          high accuracy.
    !c        On exit factr is unchanged.
    !c 
    !c      pgtol is a double precision variable.
    !c        On entry pgtol >= 0 is specified by the user.  The iteration
    !c          will stop when
    !c 
    !c                  max{|proj g_i | i = 1, ..., n} <= pgtol
    !c 
    !c          where pg_i is the ith component of the projected gradient.   
    !c        On exit pgtol is unchanged.
    !c 
    !c      wa is a double precision working array of length 
    !c        (2mmax + 5)nmax + 12mmax^2 + 12mmax.
    !c 
    !c      iwa is an integer working array of length 3nmax.
    !c 
    !c      task is a working string of characters of length 60 indicating
    !c        the current job when entering and quitting this subroutine.
    !c 
    !c      iprint is an integer variable that must be set by the user.
    !c        It controls the frequency and type of output generated:
    !c         iprint<0    no output is generated;
    !c         iprint=0    print only one line at the last iteration;
    !c         0<iprint<99 print also f and |proj g| every iprint iterations;
    !c         iprint=99   print details of every iteration except n-vectors;
    !c         iprint=100  print also the changes of active set and final x;
    !c         iprint>100  print details of every iteration including x and g;
    !c        When iprint > 0, the file iterate.dat will be created to
    !c                         summarize the iteration.
    !c 
    !c      csave is a working string of characters of length 60.
    !c 
    !c      lsave is a logical working array of dimension 4.
    !c        On exit with 'task' = NEW_X, the following information is 
    !c                                                              available:
    !c          If lsave(1) = .true.  then  the initial X has been replaced by
    !c                                      its projection in the feasible set;
    !c          If lsave(2) = .true.  then  the problem is constrained;
    !c          If lsave(3) = .true.  then  each variable has upper and lower
    !c                                      bounds;
    !c 
    !c      isave is an integer working array of dimension 44.
    !c        On exit with 'task' = NEW_X, the following information is 
    !c                                                              available:
    !c          isave(22) = the total number of intervals explored in the 
    !c                          search of Cauchy points;
    !c          isave(26) = the total number of skipped BFGS updates before 
    !c                          the current iteration;
    !c          isave(30) = the number of current iteration;
    !c          isave(31) = the total number of BFGS updates prior the current
    !c                          iteration;
    !c          isave(33) = the number of intervals explored in the search of
    !c                          Cauchy point in the current iteration;
    !c          isave(34) = the total number of function and gradient 
    !c                          evaluations;
    !c          isave(36) = the number of function value or gradient
    !c                                   evaluations in the current iteration;
    !c          if isave(37) = 0  then the subspace argmin is within the box;
    !c          if isave(37) = 1  then the subspace argmin is beyond the box;
    !c          isave(38) = the number of free variables in the current
    !c                          iteration;
    !c          isave(39) = the number of active constraints in the current
    !c                          iteration;
    !c          n + 1 - isave(40) = the number of variables leaving the set of
    !c                            active constraints in the current iteration;
    !c          isave(41) = the number of variables entering the set of active
    !c                          constraints in the current iteration.
    !c 
    !c      dsave is a double precision working array of dimension 29.
    !c        On exit with 'task' = NEW_X, the following information is
    !c                                                              available:
    !c          dsave(1) = current 'theta' in the BFGS matrix;
    !c          dsave(2) = f(x) in the previous iteration;
    !c          dsave(3) = factr*epsmch;
    !c          dsave(4) = 2-norm of the line search direction vector;
    !c          dsave(5) = the machine precision epsmch generated by the code;
    !c          dsave(7) = the accumulated time spent on searching for
    !c                                                          Cauchy points;
    !c          dsave(8) = the accumulated time spent on
    !c                                                  subspace minimization;
    !c          dsave(9) = the accumulated time spent on line search;
    !c          dsave(11) = the slope of the line search function at
    !c                                   the current point of line search;
    !c          dsave(12) = the maximum relative step length imposed in
    !c                                                            line search;
    !c          dsave(13) = the infinity norm of the projected gradient;
    !c          dsave(14) = the relative step length in the line search;
    !c          dsave(15) = the slope of the line search function at
    !c                                  the starting point of the line search;
    !c          dsave(16) = the square of the 2-norm of the line search
    !c                                                       direction vector.
    !c 
    !c      Subprograms called:
    !c 
    !c        L-BFGS-B Library ... mainlb.    
    !c 
    !c 
    !c      References:
    !c 
    !c        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
    !c        memory algorithm for bound constrained optimization'',
    !c        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
    !c 
    !c        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
    !c        limited memory FORTRAN code for solving bound constrained
    !c        optimization problems'', Tech. Report, NAM-11, EECS Department,
    !c        Northwestern University, 1994.
    !c 
    !c        (Postscript files of these papers are available via anonymous
    !c         ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
    !c 
    !c                            *  *  *
    !c 
    !c      NEOS, November 1994. (Latest revision June 1996.)
    !c      Optimization Technology Center.
    !c      Argonne National Laboratory and Northwestern University.
    !c      Written by
    !c                         Ciyou Zhu
    !c      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    !c 
    !c 
    !c      ************
    !c -jlm-jn 
    integer ::lws,lr,lz,lt,ld,lxp,lwa
    integer :: lwy,lsy,lss,lwt,lwn,lsnd

    if (task == 'start') then      
       isave = 0
       isave(1)  = m*n
       isave(2)  = m**2
       isave(3)  = 4*m**2
       isave(4)  = 1                      ! ws      m*n
       isave(5)  = isave(4)  + isave(1)   ! wy      m*n
       isave(6)  = isave(5)  + isave(1)   ! wsy     m**2
       isave(7)  = isave(6)  + isave(2)   ! wss     m**2
       isave(8)  = isave(7)  + isave(2)   ! wt      m**2
       isave(9)  = isave(8)  + isave(2)   ! wn      4*m**2
       isave(10) = isave(9)  + isave(3)   ! wsnd    4*m**2
       isave(11) = isave(10) + isave(3)   ! wz      n
       isave(12) = isave(11) + n          ! wr      n
       isave(13) = isave(12) + n          ! wd      n
       isave(14) = isave(13) + n          ! wt      n
       isave(15) = isave(14) + n          ! wxp     n
       isave(16) = isave(15) + n          ! wa      8*m    
    end if
    lws  = isave(4)
    lwy  = isave(5)
    lsy  = isave(6)
    lss  = isave(7)
    lwt  = isave(8)
    lwn  = isave(9)
    lsnd = isave(10)
    lz   = isave(11)
    lr   = isave(12)
    ld   = isave(13)
    lt   = isave(14)
    lxp  = isave(15)
    lwa  = isave(16)

    call mainlb(n,m,x,l,u,nbd,f,g,factr,pgtol,&
         wa(lws),wa(lwy),wa(lsy),wa(lss), wa(lwt), &
         wa(lwn),wa(lsnd),wa(lz),wa(lr),wa(ld),wa(lt),wa(lxp),&
         wa(lwa),&
         iwa(1),iwa(n+1),iwa(2*n+1),task,iprint, &
         csave,lsave,isave(22),dsave)           
    return
  end subroutine setulb
  !c======================= The end of setulb =============================

  subroutine mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol, ws, wy, &
       sy, ss,  wt, wn,  snd, z, r, d, t, xp, wa,&
       indx, iwhere, indx2, task, &
       iprint, csave, lsave, isave, dsave)
    implicit none
    character (len=60), intent(inout) :: task, csave
    logical, intent(inout) :: lsave(4)
    integer, intent(in) :: n, m, iprint
    integer, intent(inout)::  isave(23), indx(n), nbd(n)
    integer, intent(inout) ::  iwhere(n), indx2(n)
    real (dp), intent(in):: factr, pgtol
    real (dp), intent(inout) :: f, x(n), l(n), u(n), g(n)
    real (dp), intent(inout) :: z(n), r(n), d(n), t(n)
    !    c-jlm-jn
    real (dp), intent(inout) :: xp(n)
    real (dp), intent(inout) :: wa(8*m),  ws(n, m), wy(n, m)
    real (dp), intent(inout) :: sy(m, m), ss(m, m)
    real (dp), intent(inout) :: wt(m, m), wn(2*m, 2*m)
    real (dp), intent(inout) :: snd(2*m, 2*m), dsave(29)
    !c      ************
    !c 
    !c      Subroutine mainlb
    !c 
    !c      This subroutine solves bound constrained optimization problems by
    !c        using the compact formula of the limited memory BFGS updates.
    !c        
    !c      n is an integer variable.
    !c        On entry n is the number of variables.
    !c        On exit n is unchanged.
    !c 
    !c      m is an integer variable.
    !c        On entry m is the maximum number of variable metric
    !c           corrections allowed in the limited memory matrix.
    !c        On exit m is unchanged.
    !c 
    !c      x is a double precision array of dimension n.
    !c        On entry x is an approximation to the solution.
    !c        On exit x is the current approximation.
    !c 
    !c      l is a double precision array of dimension n.
    !c        On entry l is the lower bound of x.
    !c        On exit l is unchanged.
    !c 
    !c      u is a double precision array of dimension n.
    !c        On entry u is the upper bound of x.
    !c        On exit u is unchanged.
    !c 
    !c      nbd is an integer array of dimension n.
    !c        On entry nbd represents the type of bounds imposed on the
    !c          variables, and must be specified as follows:
    !c          nbd(i)=0 if x(i) is unbounded,
    !c                 1 if x(i) has only a lower bound,
    !c                 2 if x(i) has both lower and upper bounds,
    !c                 3 if x(i) has only an upper bound.
    !c        On exit nbd is unchanged.
    !c 
    !c      f is a double precision variable.
    !c        On first entry f is unspecified.
    !c        On final exit f is the value of the function at x.
    !c 
    !c      g is a double precision array of dimension n.
    !c        On first entry g is unspecified.
    !c        On final exit g is the value of the gradient at x.
    !c 
    !c      factr is a double precision variable.
    !c        On entry factr >= 0 is specified by the user.  The iteration
    !c          will stop when
    !c 
    !c          (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
    !c 
    !c          where epsmch is the machine precision, which is automatically
    !c          generated by the code.
    !c        On exit factr is unchanged.
    !c 
    !c      pgtol is a double precision variable.
    !c        On entry pgtol >= 0 is specified by the user.  The iteration
    !c          will stop when
    !c 
    !c                  max{|proj g_i | i = 1, ..., n} <= pgtol
    !c 
    !c          where pg_i is the ith component of the projected gradient.
    !c        On exit pgtol is unchanged.
    !c 
    !c      ws, wy, sy, and wt are double precision working arrays used to
    !c        store the following information defining the limited memory
    !c           BFGS matrix:
    !c           ws, of dimension n x m, stores S, the matrix of s-vectors;
    !c           wy, of dimension n x m, stores Y, the matrix of y-vectors;
    !c           sy, of dimension m x m, stores S'Y;
    !c           ss, of dimension m x m, stores S'S;
    !c           yy, of dimension m x m, stores Y'Y;
    !c           wt, of dimension m x m, stores the Cholesky factorization
    !c                                   of (theta*S'S+LD^(-1)L'); see eq.
    !c                                   (2.26) in [3].
    !c 
    !c      wn is a double precision working array of dimension 2m x 2m
    !c        used to store the LEL^T factorization of the indefinite matrix
    !c                  K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
    !c                      [L_a -R_z           theta*S'AA'S ]
    !c 
    !c        where     E = [-I  0]
    !c                      [ 0  I]
    !c 
    !c      snd is a double precision working array of dimension 2m x 2m
    !c        used to store the lower triangular part of
    !c                  N = [Y' ZZ'Y   L_a'+R_z']
    !c                      [L_a +R_z  S'AA'S   ]
    !c             
    !c      z(n),r(n),d(n),t(n), xp(n),wa(8*m) are double precision working arrays.
    !c        z  is used at different times to store the Cauchy point and
    !c           the Newton point.
    !c        xp is used to safeguard the projected Newton direction
    !c 
    !c      sg(m),sgo(m),yg(m),ygo(m) are double precision working arrays. 
    !c 
    !c      index is an integer working array of dimension n.
    !c        In subroutine freev, index is used to store the free and fixed
    !c           variables at the Generalized Cauchy Point (GCP).
    !c 
    !c      iwhere is an integer working array of dimension n used to record
    !c        the status of the vector x for GCP computation.
    !c        iwhere(i)=0 or -3 if x(i) is free and has bounds,
    !c                  1       if x(i) is fixed at l(i), and l(i) .ne. u(i)
    !c                  2       if x(i) is fixed at u(i), and u(i) .ne. l(i)
    !c                  3       if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
    !c                 -1       if x(i) is always free, i.e., no bounds on it.
    !c 
    !c      indx2 is an integer working array of dimension n.
    !c        Within subroutine cauchy, indx2 corresponds to the array iorder.
    !c        In subroutine freev, a list of variables entering and leaving
    !c        the free set is stored in indx2, and it is passed on to
    !c        subroutine formk with this information.
    !c 
    !c      task is a working string of characters of length 60 indicating
    !c        the current job when entering and leaving this subroutine.
    !c 
    !c      iprint is an INTEGER variable that must be set by the user.
    !c        It controls the frequency and type of output generated:
    !c         iprint<0    no output is generated;
    !c         iprint=0    print only one line at the last iteration;
    !c         0<iprint<99 print also f and |proj g| every iprint iterations;
    !c         iprint=99   print details of every iteration except n-vectors;
    !c         iprint=100  print also the changes of active set and final x;
    !c         iprint>100  print details of every iteration including x and g;
    !c        When iprint > 0, the file iterate.dat will be created to
    !c                         summarize the iteration.
    !c 
    !c      csave is a working string of characters of length 60.
    !c 
    !c      lsave is a logical working array of dimension 4.
    !c 
    !c      isave is an integer working array of dimension 23.
    !c 
    !c      dsave is a double precision working array of dimension 29.
    !c 
    !c 
    !c      Subprograms called
    !c 
    !c        L-BFGS-B Library ... cauchy, subsm, lnsrlb, formk, 
    !c 
    !c         errclb, prn1lb, prn2lb, prn3lb, active, projgr,
    !c 
    !c         freev, cmprlb, matupd, formt.
    !c 
    !c        Minpack2 Library ... timer
    !c 
    !c        Linpack Library ... dcopy, ddot.
    !c 
    !c 
    !c      References:
    !c 
    !c        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
    !c        memory algorithm for bound constrained optimization'',
    !c        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
    !c 
    !c        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
    !c        Subroutines for Large Scale Bound Constrained Optimization''
    !c        Tech. Report, NAM-11, EECS Department, Northwestern University,
    !c        1994.
    !c  
    !c        [3] R. Byrd, J. Nocedal and R. Schnabel "Representations of
    !c        Quasi-Newton Matrices and their use in Limited Memory Methods'',
    !c        Mathematical Programming 63 (1994), no. 4, pp. 129-156.
    !c 
    !c        (Postscript files of these papers are available via anonymous
    !c         ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
    !c 
    !c                            *  *  *
    !c 
    !c      NEOS, November 1994. (Latest revision June 1996.)
    !c      Optimization Technology Center.
    !c      Argonne National Laboratory and Northwestern University.
    !c      Written by
    !c                         Ciyou Zhu
    !c      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    !c 
    !c 
    !c      ************

    logical :: prjctd,cnstnd,boxed,updatd,wrk
    character (len=3) :: word
    integer :: i,k, nintol, itfile,iback,nskip
    integer :: head,col,iter,itail,iupdat
    integer :: nseg,nfgv,info,ifun
    integer :: iword,nfree,nact,ileave,nenter
    real (dp) :: theta,fold,dr,rr,tol
    real (dp) :: xstep,sbgnrm,ddum,dnorm,dtd,epsmch
    real (dp) :: cpu1,cpu2,cachyt,sbtime,lnscht,time1,time2
    real (dp) :: gd,gdold,stp,stpmx,time
    real (dp), parameter :: one=1.0d0, zero=0.0d0                

    if (task == 'start') then

       epsmch = epsilon(one)
       call timer(time1)

       !c  initialize counters and scalars when task='start'.

       !c   for the limited memory bfgs matrices:
       col    = 0  
       head   = 1   
       theta  = one  
       iupdat = 0 
       updatd = .false.
       iback  = 0  
       itail  = 0   
       iword  = 0   
       nact   = 0 
       ileave = 0
       nenter = 0

       fold   = zero 
       dnorm  = zero  
       cpu1   = zero
       gd     = zero 
       stpmx  = zero 
       sbgnrm = zero 
       stp    = zero
       gdold  = zero 
       dtd    = zero 
       xstep = zero

       !c   for operation counts:
       iter   = 0 
       nfgv   = 0 
       nseg   = 0 
       nintol = 0
       nskip  = 0 
       nfree  = n 
       ifun   = 0

       !c   for stopping tolerance:
       tol = factr*epsmch

       !c   for measuring running time:
       cachyt = 0
       sbtime = 0
       lnscht = 0

       !c  'word' records the status of subspace solutions.
       word = '---'

       !c   'info' records the termination information.
       info = 0

       ! itfile = 8
       ! if (iprint .ge. 1) then
       ! !c  open a summary file 'iterate.dat'
       ! open (8, file = 'iterate.dat', status = 'unknown')
       ! endif

       !c   check the input arguments for errors.

       call errclb(n,m,factr,l,u,nbd,task,info,k)          
       if (task(1:5) == 'error') then
          call prn3lb(n,x,f,task,iprint,info,itfile, &
               iter,nfgv,nintol,nskip,nact,sbgnrm,&
               zero,nseg,word,iback,stp,xstep,k,&
               cachyt,sbtime,lnscht)
          return
       end if

       call prn1lb(n,m,l,u,x,iprint,itfile,epsmch)

       ! c   initialize iwhere & project x onto the feasible set.

       call active(n,l,u,nbd,x,iwhere,iprint,prjctd,cnstnd,boxed)

       !c   the end of the initialization.

    else

       xstep = zero

       !c   restore local variables.   
       prjctd = lsave(1)
       cnstnd = lsave(2)
       boxed  = lsave(3)
       updatd = lsave(4)

       nintol = isave(1)
       itfile = isave(3)
       iback  = isave(4)
       nskip  = isave(5)
       head   = isave(6)
       col    = isave(7)
       itail  = isave(8)
       iter   = isave(9)
       iupdat = isave(10)
       nseg   = isave(12)
       nfgv   = isave(13)
       info   = isave(14)
       ifun   = isave(15)
       iword  = isave(16)
       nfree  = isave(17)
       nact   = isave(18)
       ileave = isave(19)
       nenter = isave(20)

       theta  = dsave(1)
       fold   = dsave(2)
       tol    = dsave(3)
       dnorm  = dsave(4)
       epsmch = dsave(5)
       cpu1   = dsave(6)
       cachyt = dsave(7)
       sbtime = dsave(8)
       lnscht = dsave(9)
       time1  = dsave(10)
       gd     = dsave(11)
       stpmx  = dsave(12)
       sbgnrm = dsave(13)
       stp    = dsave(14)
       gdold  = dsave(15)
       dtd    = dsave(16)

       !c    after returning from the driver go to the point where execution
       !c    is to resume.     

       if (task(1:5) == 'fg_ln') go to 666
       if (task(1:5) == 'new_x') go to 777
       if (task(1:5) == 'fg_st') go to 111
       if (task(1:4) == 'stop') then
          if (task(7:9) == 'cpu') then
             !c      restore the previous iterate.
             call dcopy(n,t,1,x,1)
             call dcopy(n,r,1,g,1)
             f = fold
          end if
          go to 999
       end if
    end if

    !c  compute f0 and g0.

    task = 'fg_start'
    !c   return to the driver to calculate f and g; reenter at 111.
    go to 1000

111 continue

    nfgv = 1

    !c  compute the infinity norm of the (-) projected gradient.

    call projgr(n,l,u,nbd,x,g,sbgnrm)

    if (iprint >= 1) then
       !write (6,1002) iter,f,sbgnrm
       !write (itfile,1003) iter,nfgv,sbgnrm,f
       call dblepr("at iterate: / f =  / |proj g|=", -1, (/dble(iter),f,sbgnrm/), 3)
    end if

    if (sbgnrm <= pgtol) then
       !c   terminate the algorithm.
       task = 'convergence: norm of projected gradient <= pgtol'
       go to 999
    end if

    ! ----------------- the beginning of the loop --------------------------
222 continue
    !if (iprint .ge. 99) write (6,1001) iter + 1
    if (iprint >= 99) call intpr1('iteration ', -1, iter + 1)
    iword = -1

    if (.not. cnstnd .and. col > 0) then
       !c    skip the search for gcp.
       call dcopy(n,x,1,z,1)
       wrk = updatd
       nseg = 0
       go to 333
    end if

    !cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc
    !c
    !c     compute the generalized cauchy point (gcp).
    !c
    !cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc

    call timer(cpu1)            

    call cauchy(n,x,l,u,nbd,g,indx2,iwhere,t,d,z, &
         m,wy,ws,sy,wt,theta,col,head, &
         wa(1),wa(2*m+1),wa(4*m+1),wa(6*m+1),nseg, &
         iprint,sbgnrm,info,epsmch)       

    if (info /= 0) then
       !c         singular triangular system detected; refresh the lbfgs memory.
       if(iprint >= 1) then
          call labelpr(' singular triangular system detected;', -1)
          call labelpr('   refresh the lbfgs memory and restart the iteration.', -1)
       end if
       ! if(iprint .ge. 1) write (6, 1005)
       info   = 0
       col    = 0
       head   = 1
       theta  = one
       iupdat = 0
       updatd = .false.
       call timer(cpu2)
       cachyt = cachyt + cpu2 - cpu1
       go to 222
    end if

    call timer(cpu2)
    cachyt = cachyt + cpu2 - cpu1
    nintol = nintol + nseg

    !c  count the entering and leaving variables for iter > 0; 
    !c   find the index set of free and active variables at the gcp.

    call freev(n,nfree,indx,nenter,ileave,indx2,&
         iwhere,wrk,updatd,cnstnd,iprint,iter)
    nact = n - nfree

333 continue

    !c     if there are no free variables or b=theta*i, then
    !c                             skip the subspace minimization.


    if (nfree == 0 .or. col == 0) go to 555

    !cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc
    !c
    !c     subspace minimization.
    !c
    !cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc

    call timer(cpu1)

    !c     form  the lel^t factorization of the indefinite
    !c       matrix    k = [-d -y'zz'y/theta     l_a'-r_z'  ]
    !c          [l_a -r_z           theta*s'aa's ]
    !c       where     e = [-i  0]
    !c          [ 0  i]

    if (wrk) call formk(n,nfree,indx,nenter,ileave,indx2,iupdat, &
         updatd,wn,snd,m,ws,wy,sy,theta,col,head,info)

    if (info /= 0) then
       !c  nonpositive definiteness in cholesky factorization;
       !c  refresh the lbfgs memory and restart the iteration.       
       if(iprint >= 1) then
          call labelpr(' nonpositive definiteness in cholesky factorization in formk;', -1)
          call labelpr('   refresh the lbfgs memory and restart the iteration.', -1)
       end if
       !if(iprint .ge. 1) write (6, 1006)
       info   = 0
       col    = 0
       head   = 1
       theta  = one
       iupdat = 0
       updatd = .false.
       call timer(cpu2)
       sbtime = sbtime + cpu2 - cpu1
       go to 222
    end if

    !c        compute r=-z'b(xcp-xk)-z'g (using wa(2m+1)=w'(xcp-x)
    !c                                        from 'cauchy').

    call cmprlb(n,m,x,g,ws,wy,sy,wt,z,r,wa,indx,&
         theta,col,head,nfree,cnstnd,info)

    if (info /= 0) go to 444

    !c-jlm-jn              call the direct method. 
    call subsm(n,m,nfree,indx,l,u,nbd,z,r, xp,ws,wy,&
         theta,x,g,col,head,iword,wa,wn,iprint,info)

444 continue
    if (info /= 0) then
       !c    singular triangular system detected;
       !c   refresh the lbfgs memory and restart the iteration.       
       if(iprint >= 1) then
          call labelpr(' singular triangular system detected;', -1)
          call labelpr('   refresh the lbfgs memory and restart the iteration.', -1)
       end if
       !if(iprint .ge. 1) write (6, 1005)
       info   = 0
       col    = 0
       head   = 1
       theta  = one
       iupdat = 0
       updatd = .false.
       call timer(cpu2)
       sbtime = sbtime + cpu2 - cpu1
       go to 222
    end if

    call timer(cpu2)
    sbtime = sbtime + cpu2 - cpu1

555 continue

    !cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc
    !c
    !c     line search and optimality tests.
    !c
    !cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc

    !c  generate the search direction d:=z-x.

    do i = 1, n
       d(i) = z(i) - x(i)
    end do
    call timer(cpu1)      

666 continue

    call lnsrlb(n,l,u,nbd,x,f,fold,gd,gdold,g,d,r,t,z,stp,dnorm, &
         dtd,xstep,stpmx,iter,ifun,iback,nfgv,info,task,&
         boxed,cnstnd,csave,isave(22),dsave(17))

    if (info /= 0 .or. iback >= 20) then
       ! c   restore the previous iterate.
       call dcopy(n,t,1,x,1)
       call dcopy(n,r,1,g,1)
       f = fold

       if (col == 0) then
          !c  abnormal termination.
          if (info == 0) then
             info = -9
             !c  restore the actual number of f and g evaluations etc.
             nfgv = nfgv - 1
             ifun = ifun - 1
             iback = iback - 1
          end if
          task = 'abnormal_termination_in_lnsrch'
          iter = iter + 1
          go to 999
       else
          !c  refresh the lbfgs memory and restart the iteration.
          if(iprint >= 1) then
             call labelpr(' bad direction in the line search;', -1)
             call labelpr('   refresh the lbfgs memory and restart the iteration.', -1)         
          end if
          !if(iprint .ge. 1) write (6, 1008)
          if (info == 0) nfgv = nfgv - 1
          info   = 0
          col    = 0
          head   = 1
          theta  = one
          iupdat = 0
          updatd = .false.
          task = 'restart_from_lnsrch'
          call timer(cpu2)
          lnscht = lnscht + cpu2 - cpu1
          go to 222
       end if
    else if (task(1:5) == 'fg_ln') then

       !c   return to the driver for calculating f and g; reenter at 666.
       go to 1000       
    else

       !c  calculate and print out the quantities related to the new x.
       call timer(cpu2)
       lnscht = lnscht + cpu2 - cpu1
       iter = iter + 1

       !c compute the infinity norm of the projected (-)gradient.

       call projgr(n,l,u,nbd,x,g,sbgnrm)

       !c  print iteration information.

       call prn2lb(n,x,f,g,iprint,itfile,iter,nfgv,nact,&
            sbgnrm,nseg,word,iword,iback,stp,xstep)

       go to 1000
    end if

777 continue

    !c  test for termination.

    if (sbgnrm <= pgtol) then
       !c   terminate the algorithm.
       task = 'convergence: norm of projected gradient <= pgtol'
       go to 999
    end if

    ddum = max(abs(fold), abs(f), one)
    if ((fold - f) <= tol*ddum) then
       !c   terminate the algorithm.
       task = 'convergence: rel_reduction_of_f <= factr*epsmch'
       if (iback >= 10) info = -5
       !c     i.e., to issue a warning if iback>10 in the line search.
       go to 999
    end if

    !c    compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's.
    do i = 1, n
       r(i) = g(i) - r(i)
    end do
    rr = ddot(n,r,1,r,1)

    if (stp == one) then
       dr = gd - gdold
       ddum = -gdold
    else
       dr = (gd - gdold)*stp
       call dscal(n,stp,d,1)
       ddum= -gdold*stp
    end if

    if (dr <= epsmch*ddum) then
       !c   skip the l-bfgs update.
       nskip = nskip + 1
       updatd = .false.
       if (iprint >= 1) then
          call dblepr("  nskip / ys = / -gs = / bfgs update skipped",&
               -1, (/dble(nskip), dr, ddum/), 3)
          !if (iprint .ge. 1) write (6,1004) dr, ddum
       end if
       go to 888
    end if

    !cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc
    !c
    !c     update the l-bfgs matrix.
    !c
    !cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc!cccccc

    updatd = .true.
    iupdat = iupdat + 1

    !c     update matrices ws and wy and form the middle matrix in b.
    call matupd(n,m,ws,wy,sy,ss,d,r,itail,&
         iupdat,col,head,theta,rr,dr,stp,dtd)

    !c     form the upper half of the pds t = theta*ss + l*d^(-1)*l';
    !c        store t in the upper triangular of the array wt;
    !c        cholesky factorize t to j*j' with
    !c           j' stored in the upper triangular of wt.

    call formt(m,wt,sy,ss,col,theta,info)

    if (info /= 0) then
       !c  nonpositive definiteness in cholesky factorization;
       !c   refresh the lbfgs memory and restart the iteration.
       if(iprint >= 1) then
          call labelpr(' nonpositive definiteness in cholesky factorization in formk;', -1)
          call labelpr('   refresh the lbfgs memory and restart the iteration.', -1)          
       end if
       !if(iprint .ge. 1) write (6, 1007)
       info   = 0
       col    = 0
       head   = 1
       theta  = one
       iupdat = 0
       updatd = .false.
       go to 222
    end if

    !c     now the inverse of the middle matrix in b is

    !c       [  d^(1/2)      o ] [ -d^(1/2)  d^(-1/2)*l' ]
    !c       [ -l*d^(-1/2)   j ] [  0        j'          ]

888 continue

    !c -------------------- the end of the loop -----------------------------

    goto 222

999 continue
    call timer(time2)
    time = time2 - time1
    call prn3lb(n,x,f,task,iprint,info,itfile,&
         iter,nfgv,nintol,nskip,nact,sbgnrm, &
         time,nseg,word,iback,stp,xstep,k, &
         cachyt,sbtime,lnscht)

1000 continue

    !c     save local variables.

    lsave(1) = prjctd
    lsave(2) = cnstnd
    lsave(3) = boxed
    lsave(4) = updatd

    isave(1)  = nintol
    isave(3)  = itfile
    isave(4)  = iback
    isave(5)  = nskip
    isave(6)  = head
    isave(7)  = col
    isave(8)  = itail
    isave(9)  = iter
    isave(10) = iupdat
    isave(12) = nseg
    isave(13) = nfgv
    isave(14) = info
    isave(15) = ifun
    isave(16) = iword
    isave(17) = nfree
    isave(18) = nact
    isave(19) = ileave
    isave(20) = nenter

    dsave(1)  = theta
    dsave(2)  = fold
    dsave(3)  = tol
    dsave(4)  = dnorm
    dsave(5)  = epsmch
    dsave(6)  = cpu1
    dsave(7)  = cachyt
    dsave(8)  = sbtime
    dsave(9)  = lnscht
    dsave(10) = time1
    dsave(11) = gd
    dsave(12) = stpmx
    dsave(13) = sbgnrm
    dsave(14) = stp
    dsave(15) = gdold
    dsave(16) = dtd


!!$1001 format (//,'iteration ',i5)
!!$1002 format &
!!$         (/,'at iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,d12.5)
!!$1003 format (2(1x,i4),5x,'-',5x,'-',3x,'-',5x,'-',5x,'-',8x,'-',3x,&
!!$         1p,2(1x,d10.3))
!!$1004 format ('  ys=',1p,e10.3,'  -gs=',1p,e10.3,' bfgs update skipped')
!!$1005 format (/, &
!!$         ' singular triangular system detected;',/,&
!!$         '   refresh the lbfgs memory and restart the iteration.')
!!$1006 format (/, &
!!$         ' nonpositive definiteness in cholesky factorization in formk;',/,&
!!$         '   refresh the lbfgs memory and restart the iteration.')
!!$1007 format (/, &
!!$         ' nonpositive definiteness in cholesky factorization in formt;',/,&
!!$         '   refresh the lbfgs memory and restart the iteration.')
!!$1008 format (/, &
!!$         ' bad direction in the line search;',/,&
!!$         '   refresh the lbfgs memory and restart the iteration.')

    return
  end subroutine mainlb
  !c======================= The end of mainlb =============================


  subroutine active(n, l, u, nbd, x, iwhere, iprint, prjctd, cnstnd, boxed)
    implicit none
    integer, intent(in) :: n, iprint
    real (dp), intent(inout) :: l(n), u(n)
    integer, intent(inout) :: nbd(n)
    real (dp), intent(inout) :: x(n)
    integer, intent(inout) :: iwhere(n)
    logical, intent(inout) :: prjctd, cnstnd,boxed
    !c     ************
    !c
    !c     subroutine active
    !c
    !c     this subroutine initializes iwhere and projects the initial x to
    !c       the feasible set if necessary.
    !c
    !c     iwhere is an integer array of dimension n.
    !c       on entry iwhere is unspecified.
    !c       on exit iwhere(i)=-1  if x(i) has no bounds
    !c                         3   if l(i)=u(i)
    !c                         0   otherwise.
    !c       in cauchy, iwhere is given finer gradations.
    !c
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     ************
    integer :: nbdd,i
    real (dp), parameter :: zero=0.0d0

    !c     initialize nbdd, prjctd, cnstnd and boxed.

    nbdd = 0
    prjctd = .false.
    cnstnd = .false.
    boxed = .true.

    !c     project the initial x to the easible set if necessary.

    do i = 1, n
       if (nbd(i) > 0) then
          if (nbd(i) <= 2 .and. x(i) <= l(i)) then
             if (x(i) < l(i)) then
                prjctd = .true.
                x(i) = l(i)
             end if
             nbdd = nbdd + 1
          else if (nbd(i) >= 2 .and. x(i) >= u(i)) then
             if (x(i) > u(i)) then
                prjctd = .true.
                x(i) = u(i)
             end if
             nbdd = nbdd + 1
          end if
       end if
    end do

    !c     initialize iwhere and assign values to cnstnd and boxed.

    do i = 1, n
       if (nbd(i) /= 2) boxed = .false.
       if (nbd(i) == 0) then
          !c   this variable is always free
          iwhere(i) = -1

          !c  otherwise set x(i)=mid(x(i), u(i), l(i)).
       else
          cnstnd = .true.
          if (nbd(i) == 2 .and. u(i) - l(i) <= zero) then
             !c    this variable is always fixed
             iwhere(i) = 3
          else
             iwhere(i) = 0
          end if
       end if
    end do

    if (iprint >= 0) then
       if (prjctd) then
          call labelpr( 'the initial x is infeasible. restart with its projection.', -1) 
       end if
       if (.not. cnstnd) then
          call labelpr( 'this problem is unconstrained.', -1) 
       end if
    end if

    if (iprint > 0) call intpr1('at x0 "k" variables are exactly at the bounds. k =',-1,nbdd)

!!$      if (iprint .ge. 0) then
!!$         if (prjctd) write (6,*) 'The initial X is infeasible.  Restart with its projection.'
!!$         if (.not. cnstnd) write (6,*) 'This problem is unconstrained.'
!!$      endif
!!$
!!$      if (iprint .gt. 0) write (6,1001) nbdd

!!$1001 format (/,'at x0 ',i9,' variables are exactly at the bounds')
    return
  end subroutine active
  !c======================= The end of active =============================


  subroutine bmv(m, sy, wt, col, v, p, info)
    implicit none
    integer, intent(in) :: m, col
    real (dp), intent(in) :: sy(m, m), v(2*col),  wt(m, m)
    real (dp), intent(inout) :: p(2*col)
    integer, intent(inout) :: info
    !c     ************
    !c
    !c     subroutine bmv
    !c
    !c     this subroutine computes the product of the 2m x 2m middle matrix 
    !c       in the compact l-bfgs formula of b and a 2m vector v;  
    !c       it returns the product in p.
    !c       
    !c     m is an integer variable.
    !c       on entry m is the maximum number of variable metric corrections
    !c         used to define the limited memory matrix.
    !c       on exit m is unchanged.
    !c
    !c     sy is a double precision array of dimension m x m.
    !c       on entry sy specifies the matrix s'y.
    !c       on exit sy is unchanged.
    !c
    !c     wt is a double precision array of dimension m x m.
    !c       on entry wt specifies the upper triangular matrix j' which is 
    !c         the cholesky factor of (thetas's+ld^(-1)l').
    !c       on exit wt is unchanged.
    !c
    !c     col is an integer variable.
    !c       on entry col specifies the number of s-vectors (or y-vectors)
    !c         stored in the compact l-bfgs formula.
    !c       on exit col is unchanged.
    !c
    !c     v is a double precision array of dimension 2col.
    !c       on entry v specifies vector v.
    !c       on exit v is unchanged.
    !c
    !c     p is a double precision array of dimension 2col.
    !c       on entry p is unspecified.
    !c       on exit p is the product mv.
    !c
    !c     info is an integer variable.
    !c       on entry info is unspecified.
    !c       on exit info = 0       for normal return,
    !c                    = nonzero for abnormal return when the system
    !c                                to be solved by dtrsl is singular.
    !c
    !c     subprograms called:
    !c
    !c       linpack ... dtrsl.
    !c
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     ************  
    integer :: i,k,i2
    real (dp) :: soma

    if (col == 0) return

    !c     part i: solve [  d^(1/2)      o ] [ p1 ] = [ v1 ]
    !c                   [ -l*d^(-1/2)   j ] [ p2 ]   [ v2 ].

    !c       solve jp2=v2+ld^(-1)v1.   
    p(col + 1) = v(col + 1)
    do i = 2, col
       i2 = col + i
       soma = 0.0d0
       do k = 1, i - 1
          soma = soma + sy(i,k)*v(k)/sy(k,k)
       end do
       p(i2) = v(i2) + soma
    end do
    !c     solve the triangular system
    call dtrsl(wt,m,col,p(col+1),11,info)
    if (info /= 0) return

    !c       solve d^(1/2)p1=v1.
    do i = 1, col
       p(i) = v(i)/sqrt(sy(i,i))
    end do

    !c     part ii: solve [ -d^(1/2)   d^(-1/2)*l'  ] [ p1 ] = [ p1 ]
    !c                    [  0         j'           ] [ p2 ]   [ p2 ]. 

    !c       solve j^tp2=p2.
    call dtrsl(wt,m,col,p(col+1),1,info)
    if (info /= 0) return

    !c       compute p1=-d^(-1/2)(p1-d^(-1/2)l'p2)
    !c                 =-d^(-1/2)p1+d^(-1)l'p2.    
    do i = 1, col
       p(i) = -p(i)/sqrt(sy(i,i))
    end do
    do i = 1, col
       soma = 0.d0
       do k = i + 1, col
          soma = soma + sy(k,i)*p(col+k)/sy(i,i)
       end do
       p(i) = p(i) + soma
    end do
    return
  end subroutine bmv
  !c======================== The end of bmv ===============================

  subroutine cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp,&
       m, wy, ws, sy, wt,  theta, col, head, p, c, wbp,&
       v, nseg,  iprint, sbgnrm, info, epsmch)
    implicit none
    integer, intent(in) :: n, m, iprint, nbd(n)
    integer, intent(in) :: col, head
    integer, intent(inout) :: iorder(n), nseg
    integer, intent(inout) :: iwhere(n), info
    real (dp), intent(in) :: x(n), l(n), u(n), g(n), sbgnrm, epsmch
    real (dp), intent(in) :: ws(n, col), wy(n, col),  sy(m, m), wt(m, m), theta
    real (dp), intent(inout) :: t(n), d(n), xcp(n), p(2*m), c(2*m), wbp(2*m),  v(2*m)
    !c     ************
    !c
    !c     subroutine cauchy
    !c
    !c     for given x, l, u, g (with sbgnrm > 0), and a limited memory
    !c       bfgs matrix b defined in terms of matrices wy, ws, wt, and
    !c       scalars head, col, and theta, this subroutine computes the
    !c       generalized cauchy point (gcp), defined as the first local
    !c       minimizer of the quadratic
    !c
    !c                  q(x + s) = g's + 1/2 s'bs
    !c
    !c       along the projected gradient direction p(x-tg,l,u).
    !c       the routine returns the gcp in xcp. 
    !c       
    !c     n is an integer variable.
    !c       on entry n is the dimension of the problem.
    !c       on exit n is unchanged.
    !c
    !c     x is a double precision array of dimension n.
    !c       on entry x is the starting point for the gcp computation.
    !c       on exit x is unchanged.
    !c
    !c     l is a double precision array of dimension n.
    !c       on entry l is the lower bound of x.
    !c       on exit l is unchanged.
    !c
    !c     u is a double precision array of dimension n.
    !c       on entry u is the upper bound of x.
    !c       on exit u is unchanged.
    !c
    !c     nbd is an integer array of dimension n.
    !c       on entry nbd represents the type of bounds imposed on the
    !c         variables, and must be specified as follows:
    !c         nbd(i)=0 if x(i) is unbounded,
    !c                1 if x(i) has only a lower bound,
    !c                2 if x(i) has both lower and upper bounds, and
    !c                3 if x(i) has only an upper bound. 
    !c       on exit nbd is unchanged.
    !c
    !c     g is a double precision array of dimension n.
    !c       on entry g is the gradient of f(x).  g must be a nonzero vector.
    !c       on exit g is unchanged.
    !c
    !c     iorder is an integer working array of dimension n.
    !c       iorder will be used to store the breakpoints in the piecewise
    !c       linear path and free variables encountered. on exit,
    !c         iorder(1),...,iorder(nleft) are indices of breakpoints
    !c                                which have not been encountered; 
    !c         iorder(nleft+1),...,iorder(nbreak) are indices of
    !c                                     encountered breakpoints; and
    !c         iorder(nfree),...,iorder(n) are indices of variables which
    !c                 have no bound constraits along the search direction.
    !c
    !c     iwhere is an integer array of dimension n.
    !c       on entry iwhere indicates only the permanently fixed (iwhere=3)
    !c       or free (iwhere= -1) components of x.
    !c       on exit iwhere records the status of the current x variables.
    !c       iwhere(i)=-3  if x(i) is free and has bounds, but is not moved
    !c                 0   if x(i) is free and has bounds, and is moved
    !c                 1   if x(i) is fixed at l(i), and l(i) .ne. u(i)
    !c                 2   if x(i) is fixed at u(i), and u(i) .ne. l(i)
    !c                 3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
    !c                 -1  if x(i) is always free, i.e., it has no bounds.
    !c
    !c     t is a double precision working array of dimension n. 
    !c       t will be used to store the break points.
    !c
    !c     d is a double precision array of dimension n used to store
    !c       the cauchy direction p(x-tg)-x.
    !c
    !c     xcp is a double precision array of dimension n used to return the
    !c       gcp on exit.
    !c
    !c     m is an integer variable.
    !c       on entry m is the maximum number of variable metric corrections 
    !c         used to define the limited memory matrix.
    !c       on exit m is unchanged.
    !c
    !c     ws, wy, sy, and wt are double precision arrays.
    !c       on entry they store information that defines the
    !c                             limited memory bfgs matrix:
    !c         ws(n,m) stores s, a set of s-vectors;
    !c         wy(n,m) stores y, a set of y-vectors;
    !c         sy(m,m) stores s'y;
    !c         wt(m,m) stores the
    !c                 cholesky factorization of (theta*s's+ld^(-1)l').
    !c       on exit these arrays are unchanged.
    !c
    !c     theta is a double precision variable.
    !c       on entry theta is the scaling factor specifying b_0 = theta i.
    !c       on exit theta is unchanged.
    !c
    !c     col is an integer variable.
    !c       on entry col is the actual number of variable metric
    !c         corrections stored so far.
    !c       on exit col is unchanged.
    !c
    !c     head is an integer variable.
    !c       on entry head is the location of the first s-vector (or y-vector)
    !c         in s (or y).
    !c       on exit col is unchanged.
    !c
    !c     p is a double precision working array of dimension 2m.
    !c       p will be used to store the vector p = w^(t)d.
    !c
    !c     !c is a double precision working array of dimension 2m.
    !c       !c will be used to store the vector !c = w^(t)(xcp-x).
    !c
    !c     wbp is a double precision working array of dimension 2m.
    !c       wbp will be used to store the row of w corresponding
    !c         to a breakpoint.
    !c
    !c     v is a double precision working array of dimension 2m.
    !c
    !c     nseg is an integer variable.
    !c       on exit nseg records the number of quadrati!c segments explored
    !c         in searching for the gcp.
    !c
    !c     sg and yg are double precision arrays of dimension m.
    !c       on entry sg  and yg store s'g and y'g correspondingly.
    !c       on exit they are unchanged. 
    !c 
    !c     iprint is an integer variable that must be set by the user.
    !c       it controls the frequency and type of output generated:
    !c        iprint<0    no output is generated;
    !c        iprint=0    print only one line at the last iteration;
    !c        0<iprint<99 print also f and |proj g| every iprint iterations;
    !c        iprint=99   print details of every iteration except n-vectors;
    !c        iprint=100  print also the changes of active set and final x;
    !c        iprint>100  print details of every iteration including x and g;
    !c       when iprint > 0, the file iterate.dat will be created to
    !c                        summarize the iteration.
    !c
    !c     sbgnrm is a double precision variable.
    !c       on entry sbgnrm is the norm of the projected gradient at x.
    !c       on exit sbgnrm is unchanged.
    !c
    !c     info is an integer variable.
    !c       on entry info is 0.
    !c       on exit info = 0       for normal return,
    !c                    = nonzero for abnormal return when the the system
    !c                              used in routine bmv is singular.
    !c
    !c     subprograms called:
    !c 
    !c       l-bfgs-b library ... hpsolb, bmv.
    !c
    !c       linpack ... dscal dcopy, daxpy.
    !c
    !c
    !c     references:
    !c
    !c       [1] r. h. byrd, p. lu, j. nocedal and c. zhu, ``a limited
    !c       memory algorithm for bound constrained optimization'',
    !c       siam j. scientifi!!!c     omputing 16 (1995), no. 5, pp. 1190--1208.
    !c
    !c       [2] c. zhu, r.h. byrd, p. lu, j. nocedal, ``l-bfgs-b: fortran
    !c       subroutines for large scale bound constrained optimization''
    !c       tech. report, nam-11, eecs department, northwestern university,
    !c       1994.
    !c
    !c       (postscript files of these papers are available via anonymous
    !c        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     ************
    logical :: xlower,xupper,bnded
    integer :: i,j,col2,nfree,nbreak,pointr, ibp,nleft,ibkmin,iter
    real (dp) :: f1,f2,dt,dtm,tsum,dibp,zibp,dibp2,bkmin
    real (dp) :: tu,tl,wmc,wmp,wmw,tj,tj0,neggi, f2_org
    real (dp), parameter :: one=1.0d0, zero=0.0d0

    !c     check the status of the variables, reset iwhere(i) if necessary;
    !c       compute the cauchy direction d and the breakpoints t; initialize
    !c       the derivative f1 and the vector p = w'd (for theta = 1).

    if (sbgnrm <= zero) then
       if (iprint >= 0) call labelpr('subgnorm = 0.  gcp = x.', -1)
       !if (iprint .ge. 0) write (6,*) 'Subgnorm = 0.  GCP = X.'
       call dcopy(n,x,1,xcp,1)
       return
    end if

    bnded = .true.
    nfree = n + 1
    nbreak = 0
    ibkmin = 0
    bkmin = zero
    col2 = 2*col
    f1= zero
    if (iprint >= 99) call labelpr('---------------- cauchy entered-------------------', -1)
    !if (iprint .ge. 99) write (6,3010)

    !c     we set p to zero and build it up as we determine d.
    do i = 1, col2
       p(i) = zero
    end do

    !c     in the following loop we determine for each variable its bound
    !c        status and its breakpoint, and update p accordingly.
    !c        smallest breakpoint is identified.

    do i = 1, n
       neggi = -g(i)
       if (iwhere(i) /= 3 .and. iwhere(i) /= -1) then          
          !c             if x(i) is not a constant and has bounds,
          !c             compute the difference between x(i) and its bounds.       
          if (nbd(i) <= 2) tl = x(i) - l(i)
          if (nbd(i) >= 2) tu = u(i) - x(i)

          !c           if a variable is close enough to a bound
          !c             we treat it as at bound.          
          xlower = nbd(i) <= 2 .and. tl <= zero
          xupper = nbd(i) >= 2 .and. tu <= zero

          !c              reset iwhere(i).
          iwhere(i) = 0
          if (xlower) then
             if (neggi <= zero) iwhere(i) = 1
          else if (xupper) then
             if (neggi >= zero) iwhere(i) = 2
          else
             if (abs(neggi) <= zero) iwhere(i) = -3
          end if
       end if

       pointr = head
       if (iwhere(i) /= 0 .and. iwhere(i) /= -1) then
          d(i) = zero
       else
          d(i) = neggi
          f1 = f1 - neggi*neggi
          !c             calculate p := p - w'e_i* (g_i).
          do j = 1, col
             p(j) = p(j) + wy(i,pointr)* neggi
             p(col + j) = p(col + j) + ws(i,pointr)*neggi
             pointr = mod(pointr,m) + 1
          end do
          if (nbd(i) <= 2 .and. nbd(i) /= 0 .and. neggi < zero) then
             !c     x(i) + d(i) is bounded; compute t(i).
             nbreak = nbreak + 1
             iorder(nbreak) = i
             t(nbreak) = tl/(-neggi)
             if (nbreak == 1 .or. t(nbreak) < bkmin) then
                bkmin = t(nbreak)
                ibkmin = nbreak
             end if
          else if (nbd(i) >= 2 .and. neggi > zero) then
             !c         x(i) + d(i) is bounded; compute t(i).
             nbreak = nbreak + 1
             iorder(nbreak) = i
             t(nbreak) = tu/neggi
             if (nbreak == 1 .or. t(nbreak) < bkmin) then
                bkmin = t(nbreak)
                ibkmin = nbreak
             end if
          else
             !c       x(i) + d(i) is not bounded.
             nfree = nfree - 1
             iorder(nfree) = i
             if (abs(neggi) > zero) bnded = .false.
          end if
       end if
    end do

    !c     the indices of the nonzero components of d are now stored
    !c       in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n).
    !c       the smallest of the nbreak breakpoints is in t(ibkmin)=bkmin.

    if (theta /= one) then
       !c                   complete the initialization of p for theta not= one.
       call dscal(col,theta,p(col+1),1)
    end if

    !c     initialize gcp xcp = x.

    call dcopy(n,x,1,xcp,1)

    if (nbreak == 0 .and. nfree == n + 1) then
       !c                  is a zero vector, return with the initial xcp as gcp.
       return
    end if

    !c     initialize c = w'(xcp - x) = 0.

    do j = 1, col2
       c(j) = zero
    end do

    !c     initialize derivative f2.

    f2 = -theta*f1
    f2_org = f2
    if (col > 0) then
       call bmv(m,sy,wt,col,p,v,info)
       if (info /= 0) return
       f2 = f2 - ddot(col2,v,1,p,1)
    end if
    dtm = -f1/f2
    tsum = zero
    nseg = 1
    if (iprint >= 99) call intpr1('there are nbreak breakpoints. nbreak = ', -1, nbreak)
    !      if (iprint .ge. 99) write (6,*) 'There are ',nbreak,'  breakpoints '

    !c     if there are no breakpoints, locate the gcp and return.

    if (nbreak == 0) go to 888

    nleft = nbreak
    iter = 1
    tj = zero

    !------------------- the beginning of the loop -------------------------

777 continue

    !c     find the next smallest breakpoint;
    !c       compute dt = t(nleft) - t(nleft + 1).

    tj0 = tj
    if (iter == 1) then
       !c         since we already have the smallest breakpoint we need not do
       !c         heapsort yet. often only one breakpoint is used and the
       !c         cost of heapsort is avoided.
       tj = bkmin
       ibp = iorder(ibkmin)
    else
       if (iter == 2) then
          !c             replace the already used smallest breakpoint with the
          !c             breakpoint numbered nbreak > nlast, before heapsort call.          
          if (ibkmin /= nbreak) then
             t(ibkmin) = t(nbreak)
             iorder(ibkmin) = iorder(nbreak)
          end if
          !c        update heap structure of breakpoints
          !c           (if iter=2, initialize heap).          
       end if
       call hpsolb(nleft,t,iorder,iter-2)
       tj = t(nleft)
       ibp = iorder(nleft)
    end if

    dt = tj - tj0

    ! if (dt /= zero .and. iprint >= 100) then
    ! write (*,4011) nseg,f1,f2
    ! write (*,5010) dt
    ! write (*,6010) dtm
    ! end if

    !c     if a minimizer is within this interval, locate the gcp and return.

    if (dtm < dt) go to 888

    !c     otherwise fix one variable and
    !c       reset the corresponding component of d to zero.

    tsum = tsum + dt
    nleft = nleft - 1
    iter = iter + 1
    dibp = d(ibp)
    d(ibp) = zero

    if (dibp > zero) then
       zibp = u(ibp) - x(ibp)
       xcp(ibp) = u(ibp)
       iwhere(ibp) = 2
    else
       zibp = l(ibp) - x(ibp)
       xcp(ibp) = l(ibp)
       iwhere(ibp) = 1
    end if

    ! if (iprint >= 100) write (*,*) 'variable ',ibp,' is fixed.'
    if (nleft == 0 .and. nbreak == n) then
       !c    all n variables are fixed,
       !c   return with xcp as gcp.       
       dtm = dt
       go to 999
    end if

    !c     update the derivative information.

    nseg = nseg + 1
    dibp2 = dibp**2

    !c     update f1 and f2.

    !c        temporarily set f1 and f2 for col=0.    
    f1 = f1 + dt*f2 + dibp2 - theta*dibp*zibp
    f2 = f2 - theta*dibp2

    if (col > 0) then
       !c    update c = c + dt*p.
       call daxpy(col2,dt,p,1,c,1)

       ! c   choose wbp,
       !c    the row of w corresponding to the breakpoint encountered.
       pointr = head
       do j = 1,col
          wbp(j) = wy(ibp,pointr)
          wbp(col + j) = theta*ws(ibp,pointr)
          pointr = mod(pointr,m) + 1
       end do

       !c    compute (wbp)mc, (wbp)mp, and (wbp)m(wbp)'.
       call bmv(m,sy,wt,col,wbp,v,info)
       if (info /= 0) return
       wmc = ddot(col2,c,1,v,1)
       wmp = ddot(col2,p,1,v,1)
       wmw = ddot(col2,wbp,1,v,1)

       !c           update p = p - dibp*wbp. 
       call daxpy(col2,-dibp,wbp,1,p,1)

       !c           complete updating f1 and f2 while col > 0.
       f1 = f1 + dibp*wmc
       f2 = f2 + 2.0d0*dibp*wmp -dibp2*wmw
    end if

    f2 = max(epsmch*f2_org,f2)
    if (nleft > 0) then
       dtm = -f1/f2
       go to 777
       !c                 to repeat the loop for unsearched intervals
    else if(bnded) then
       f1 = zero
       f2 = zero
       dtm = zero
    else
       dtm = -f1/f2
    end if

    !------------------- the end of the loop -------------------------------

888 continue


    ! if (iprint >= 99) then
    ! write (*,*)
    ! write (*,*) 'gcp found in this segment'
    ! write (*,4010) nseg,f1,f2
    ! write (*,6010) dtm
    ! end if
    if (dtm <= zero) dtm = zero
    tsum = tsum + dtm

    !c     move free variables (i.e., the ones w/o breakpoints) and 
    !c       the variables whose breakpoints haven't been reached
    call daxpy(n,tsum,d,1,xcp,1)

999 continue

    !c     update c = c + dtm*p = w'(x^c - x) 
    !c       which will be used in computing r = z'(b(x^c - x) + g).

    if (col > 0) call daxpy(col2,dtm,p,1,c,1)
    ! if (iprint > 100) write (*,1010) (xcp(i),i = 1,n)
    ! if (iprint >= 99) write (*,2010)

!!$1010 format ('cauchy x = ',/,(4x,1p,6(1x,d11.4)))
!!$2010 format (/,'---------------- exit cauchy----------------------',/)
!!$3010 format (/,'---------------- cauchy entered-------------------')
!!$4010 format ('piece ',i3,' --f1, f2 at start point ',1p,2(1x,d11.4))
!!$4011 format (/,'piece ',i3,' --f1, f2 at start point ', 1p,2(1x,d11.4))
!!$5010 format ('distance to the next break point = ',1p,d11.4)
!!$6010 format ('distance to the stationary point = ',1p,d11.4)
    return
  end subroutine cauchy
  !c====================== The end of cauchy ==============================

  subroutine cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, indx, &
       theta, col, head, nfree, cnstnd, info)
    implicit none
    integer, intent(in) :: n, m, col, head, nfree
    real (dp), intent(in) :: x(n), g(n), wy(n, m),ws(n, m), theta
    real (dp), intent(inout) :: sy(m, m), wt(m, m), z(n), wa(4*m)
    real (dp), intent(inout) :: r(n)
    integer, intent(in) :: indx(n)
    logical, intent(in) :: cnstnd
    integer, intent(inout) :: info
    !c     ************
    !c
    !c     subroutine cmprlb 
    !c
    !c       this subroutine computes r=-z'b(xcp-xk)-z'g by using 
    !c         wa(2m+1)=w'(xcp-x) from subroutine cauchy.
    !c
    !c     subprograms called:
    !c
    !c       l-bfgs-b library ... bmv.
    !c
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     ************
    integer :: i,j,k,pointr
    real (dp) :: a1,a2

    if (.not. cnstnd .and. col > 0) then
       do i = 1, n
          r(i) = -g(i)
       end do
    else
       do i = 1, nfree
          k = indx(i)
          r(i) = -theta*(z(k) - x(k)) - g(k)
       end do
       call bmv(m,sy,wt,col,wa(2*m+1),wa(1),info)
       if (info /= 0) then
          info = -8
          return
       end if
       pointr = head
       do j = 1, col
          a1 = wa(j)
          a2 = theta*wa(col + j)
          do i = 1, nfree
             k = indx(i)
             r(i) = r(i) + wy(k,pointr)*a1 + ws(k,pointr)*a2
          end do
          pointr = mod(pointr,m) + 1
       end do
    end if
    return
  end subroutine cmprlb
  !c======================= The end of cmprlb =============================

  subroutine errclb(n, m, factr, l, u, nbd, task, info, k)
    implicit none
    integer, intent(in) :: n, m, nbd(n)
    real (dp), intent(in) :: factr, u(n), l(n)
    character (len=60), intent(inout):: task
    integer, intent(inout) :: info, k
    !c     ************
    !c
    !c     subroutine errclb
    !c
    !c     this subroutine checks the validity of the input data.
    !c
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     ************    
    integer :: i
    real (dp), parameter :: one=1.0d0, zero=0.0d0


    !c     check the input arguments for errors.

    if (n <= 0) task = 'error: n .le. 0'
    if (m <= 0) task = 'error: m .le. 0'
    if (factr < zero) task = 'error: factr .lt. 0'

    !c     check the validity of the arrays nbd(i), u(i), and l(i).
    do i = 1, n
       if (nbd(i) < 0 .or. nbd(i) > 3) then
          task = 'error: invalid nbd'
          info = -6
          k = i
       end if
       if (nbd(i) == 2) then
          if (l(i) > u(i)) then
             task = 'error: no feasible solution'
             info = -7
             k = i
          end if
       end if
    end do
    return
  end subroutine errclb
  !c======================= The end of errclb =============================

  subroutine formk(n, nsub, ind, nenter, ileave, indx2, iupdat,&
       updatd, wn, wn1, m, ws, wy, sy, theta, col,&
       head, info)    
    implicit none
    integer, intent(in) :: n, m, nsub, ind(n), nenter, ileave, indx2(n), col, head
    integer, intent(inout) :: iupdat
    logical, intent(inout) :: updatd
    real (dp), intent(inout) :: wn(2*m, 2*m), wn1(2*m, 2*m)
    real (dp), intent(in) :: ws(n, m), wy(n, m), sy(m, m), theta
    integer, intent(inout) :: info
    !c     ************
    !c
    !c     subroutine formk 
    !c
    !c     this subroutine forms  the lel^t factorization of the indefinite
    !c
    !c       matrix    k = [-d -y'zz'y/theta     l_a'-r_z'  ]
    !c                     [l_a -r_z           theta*s'aa's ]
    !c                                                    where e = [-i  0]
    !c                                                              [ 0  i]
    !c     the matrix k can be shown to be equal to the matrix m^[-1]n
    !c       occurring in section 5.1 of [1], as well as to the matrix
    !c       mbar^[-1] nbar in section 5.3.
    !c
    !c     n is an integer variable.
    !c       on entry n is the dimension of the problem.
    !c       on exit n is unchanged.
    !c
    !c     nsub is an integer variable
    !c       on entry nsub is the number of subspace variables in free set.
    !c       on exit nsub is not changed.
    !c
    !c     ind is an integer array of dimension nsub.
    !c       on entry ind specifies the indices of subspace variables.
    !c       on exit ind is unchanged. 
    !c
    !c     nenter is an integer variable.
    !c       on entry nenter is the number of variables entering the 
    !c         free set.
    !c       on exit nenter is unchanged. 
    !c
    !c     ileave is an integer variable.
    !c       on entry indx2(ileave),...,indx2(n) are the variables leaving
    !c         the free set.
    !c       on exit ileave is unchanged. 
    !c
    !c     indx2 is an integer array of dimension n.
    !c       on entry indx2(1),...,indx2(nenter) are the variables entering
    !c         the free set, while indx2(ileave),...,indx2(n) are the
    !c         variables leaving the free set.
    !c       on exit indx2 is unchanged. 
    !c
    !c     iupdat is an integer variable.
    !c       on entry iupdat is the total number of bfgs updates made so far.
    !c       on exit iupdat is unchanged. 
    !c
    !c     updatd is a logical variable.
    !c       on entry 'updatd' is true if the l-bfgs matrix is updatd.
    !c       on exit 'updatd' is unchanged. 
    !c
    !c     wn is a double precision array of dimension 2m x 2m.
    !c       on entry wn is unspecified.
    !c       on exit the upper triangle of wn stores the lel^t factorization
    !c         of the 2*col x 2*col indefinite matrix
    !c                     [-d -y'zz'y/theta     l_a'-r_z'  ]
    !c                     [l_a -r_z           theta*s'aa's ]
    !c
    !c     wn1 is a double precision array of dimension 2m x 2m.
    !c       on entry wn1 stores the lower triangular part of 
    !c                     [y' zz'y   l_a'+r_z']
    !c                     [l_a+r_z   s'aa's   ]
    !c         in the previous iteration.
    !c       on exit wn1 stores the corresponding updated matrices.
    !c       the purpose of wn1 is just to store these inner products
    !c       so they can be easily updated and inserted into wn.
    !c
    !c     m is an integer variable.
    !c       on entry m is the maximum number of variable metric corrections
    !c         used to define the limited memory matrix.
    !c       on exit m is unchanged.
    !c
    !c     ws, wy, sy, and wtyy are double precision arrays;
    !c     theta is a double precision variable;
    !c     col is an integer variable;
    !c     head is an integer variable.
    !c       on entry they store the information defining the
    !c                                          limited memory bfgs matrix:
    !c         ws(n,m) stores s, a set of s-vectors;
    !c         wy(n,m) stores y, a set of y-vectors;
    !c         sy(m,m) stores s'y;
    !c         wtyy(m,m) stores the cholesky factorization
    !c                                   of (theta*s's+ld^(-1)l')
    !c         theta is the scaling factor specifying b_0 = theta i;
    !c         col is the number of variable metric corrections stored;
    !c         head is the location of the 1st s- (or y-) vector in s (or y).
    !c       on exit they are unchanged.
    !c
    !c     info is an integer variable.
    !c       on entry info is unspecified.
    !c       on exit info =  0 for normal return;
    !c                    = -1 when the 1st cholesky factorization failed;
    !c                    = -2 when the 2st cholesky factorization failed.
    !c
    !c     subprograms called:
    !c
    !c       linpack ... dcopy, dpofa, dtrsl.
    !c
    !c
    !c     references:
    !c       [1] r. h. byrd, p. lu, j. nocedal and c. zhu, ``a limited
    !c       memory algorithm for bound constrained optimization'',
    !c       siam j. scientifi!!!c     omputing 16 (1995), no. 5, pp. 1190--1208.
    !c
    !c       [2] c. zhu, r.h. byrd, p. lu, j. nocedal, ``l-bfgs-b: a
    !c       limited memory fortran code for solving bound constrained
    !c       optimization problems'', tech. report, nam-11, eecs department,
    !c       northwestern university, 1994.
    !c
    !c       (postscript files of these papers are available via anonymous
    !c        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     ************    
    integer :: m2,ipntr,jpntr,iy,is,jy,js,is1,js1,k1,i,k
    integer :: col2,pbegin,pend,dbegin,dend,upcl
    real (dp) :: temp1,temp2,temp3,temp4
    real (dp), parameter :: one=1.0d0, zero=0.0d0

    !c     form the lower triangular part of
    !c               wn1 = [y' zz'y   l_a'+r_z'] 
    !c                     [l_a+r_z   s'aa's   ]
    !c        where l_a is the strictly lower triangular part of s'aa'y
    !c              r_z is the upper triangular part of s'zz'y.

    if (updatd) then
       if (iupdat > m) then
          !c                shift old part of wn1.
          do jy = 1, m - 1
             js = m + jy
             call dcopy(m-jy,wn1(jy+1,jy+1),1,wn1(jy,jy),1)
             call dcopy(m-jy,wn1(js+1,js+1),1,wn1(js,js),1)
             call dcopy(m-1,wn1(m+2,jy+1),1,wn1(m+1,jy),1)
          end do
       end if

       !c          put new rows in blocks (1,1), (2,1) and (2,2).
       pbegin = 1
       pend = nsub
       dbegin = nsub + 1
       dend = n
       iy = col
       is = m + col
       ipntr = head + col - 1
       if (ipntr > m) ipntr = ipntr - m
       jpntr = head
       do jy = 1, col
          js = m + jy
          temp1 = zero
          temp2 = zero
          temp3 = zero
          !c             compute element jy of row 'col' of y'zz'y
          do k = pbegin, pend
             k1 = ind(k)
             temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
          end do
          !c             compute elements jy of row 'col' of l_a and s'aa's
          do k = dbegin, dend
             k1 = ind(k)
             temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
             temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
          end do
          wn1(iy,jy) = temp1
          wn1(is,js) = temp2
          wn1(is,jy) = temp3
          jpntr = mod(jpntr,m) + 1
       end do

       !c          put new column in block (2,1).
       jy = col
       jpntr = head + col - 1
       if (jpntr > m) jpntr = jpntr - m
       ipntr = head
       do i = 1, col
          is = m + i
          temp3 = zero
          !c             compute element i of column 'col' of r_z
          do k = pbegin, pend
             k1 = ind(k)
             temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
          end do
          ipntr = mod(ipntr,m) + 1
          wn1(is,jy) = temp3
       end do
       upcl = col - 1
    else
       upcl = col
    end if

    !c       modify the old parts in blocks (1,1) and (2,2) due to changes
    !c       in the set of free variables.
    ipntr = head
    do iy = 1, upcl
       is = m + iy
       jpntr = head
       do jy = 1, iy
          js = m + jy
          temp1 = zero
          temp2 = zero
          temp3 = zero
          temp4 = zero
          do k = 1, nenter
             k1 = indx2(k)
             temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
             temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
          end do
          do k = ileave, n
             k1 = indx2(k)
             temp3 = temp3 + wy(k1,ipntr)*wy(k1,jpntr)
             temp4 = temp4 + ws(k1,ipntr)*ws(k1,jpntr)
          end do
          wn1(iy,jy) = wn1(iy,jy) + temp1 - temp3
          wn1(is,js) = wn1(is,js) - temp2 + temp4
          jpntr = mod(jpntr,m) + 1
       end do
       ipntr = mod(ipntr,m) + 1
    end do

    !c       modify the old parts in block (2,1).
    ipntr = head
    do is = m + 1, m + upcl
       jpntr = head
       do jy = 1, upcl
          temp1 = zero
          temp3 = zero
          do k = 1, nenter
             k1 = indx2(k)
             temp1 = temp1 + ws(k1,ipntr)*wy(k1,jpntr)
          end do
          do k = ileave, n
             k1 = indx2(k)
             temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
          end do
          if (is <= jy + m) then
             wn1(is,jy) = wn1(is,jy) + temp1 - temp3
          else
             wn1(is,jy) = wn1(is,jy) - temp1 + temp3
          end if
          jpntr = mod(jpntr,m) + 1
       end do
       ipntr = mod(ipntr,m) + 1
    end do

    !c     form the upper triangle of wn = [d+y' zz'y/theta   -l_a'+r_z' ] 
    !c                                     [-l_a +r_z        s'aa's*theta]

    m2 = 2*m
    do iy = 1, col
       is = col + iy
       is1 = m + iy
       do jy = 1, iy
          js = col + jy
          js1 = m + jy
          wn(jy,iy) = wn1(iy,jy)/theta
          wn(js,is) = wn1(is1,js1)*theta
       end do
       do jy = 1, iy - 1
          wn(jy,is) = -wn1(is1,jy)
       end do
       do jy = iy, col
          wn(jy,is) = wn1(is1,jy)
       end do
       wn(iy,iy) = wn(iy,iy) + sy(iy,iy)
    end do

    !c     form the upper triangle of wn= [  ll'            l^-1(-l_a'+r_z')] 
    !c                                    [(-l_a +r_z)l'^-1   s'aa's*theta  ]

    !c        first cholesky factor (1,1) block of wn to get ll'
    !c                          with l' stored in the upper triangle of wn.
    call dpofa(wn,m2,col,info)
    if (info /= 0) then
       info = -1
       return
    end if
    !c        then form l^-1(-l_a'+r_z') in the (1,2) block.
    col2 = 2*col
    do js = col+1 ,col2
       call dtrsl(wn,m2,col,wn(1,js),11,info)
    end do

    !c     form s'aa's*theta + (l^-1(-l_a'+r_z'))'l^-1(-l_a'+r_z') in the
    !c        upper triangle of (2,2) block of wn.

    do is = col+1, col2
       do js = is, col2
          wn(is,js) = wn(is,js) + ddot(col,wn(1,is),1,wn(1,js),1)
       end do
    end do

    !c     cholesky factorization of (2,2) block of wn.

    call dpofa(wn(col+1,col+1),m2,col,info)
    if (info /= 0) then
       info = -2
       return
    end if
    return
  end subroutine formk
  !c======================= The end of formk ==============================

  subroutine formt(m, wt, sy, ss, col, theta, info)
    implicit none
    integer, intent(in) :: m, col
    real (dp), intent(inout) :: wt(m, m)
    real (dp), intent(in) :: sy(m, m), ss(m, m), theta
    integer, intent(inout) :: info
    !c     ************
    !c
    !c     subroutine formt
    !c
    !c       this subroutine forms the upper half of the pos. def. and symm.
    !c         t = theta*ss + l*d^(-1)*l', stores t in the upper triangle
    !c         of the array wt, and performs the cholesky factorization of t
    !c         to produce j*j', with j' stored in the upper triangle of wt.
    !c
    !c     subprograms called:
    !c
    !c       linpack ... dpofa.
    !c
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     ************
    integer :: i,j,k,k1
    real (dp) :: ddum
    real (dp), parameter :: zero=0.0d0

    !c     form the upper half of  t = theta*ss + l*d^(-1)*l',
    !c        store t in the upper triangle of the array wt.

    do j = 1, col
       wt(1,j) = theta*ss(1,j)
    end do
    do i = 2, col
       do j = i, col
          k1 = min(i,j) - 1
          ddum = zero
          do k = 1, k1
             ddum = ddum + sy(i,k)*sy(j,k)/sy(k,k)
          end do
          wt(i,j) = ddum + theta*ss(i,j)
       end do
    end do

    !c     cholesky factorize t to j*j' with 
    !c        j' stored in the upper triangle of wt.

    call dpofa(wt,m,col,info)
    if (info /= 0) then
       info = -3
    end if
    return
  end subroutine formt
  !c======================= The end of formt ==============================

  subroutine freev(n, nfree, indx, nenter, ileave, indx2, &
       iwhere, wrk, updatd, cnstnd,iprint, iter)
    implicit none
    integer, intent(in) :: n, iprint
    integer, intent(inout) :: nfree, indx(n), iwhere(n), iter
    integer, intent(inout) :: nenter, ileave
    integer, intent(inout) :: indx2(n)
    logical, intent(inout) :: wrk
    logical, intent(in) :: updatd
    logical, intent(inout) :: cnstnd
    !c     ************
    !c
    !c     subroutine freev 
    !c
    !c     this subroutine counts the entering and leaving variables when
    !c       iter > 0, and finds the index set of free and active variables
    !c       at the gcp.
    !c
    !c     cnstnd is a logical variable indicating whether bounds are present
    !c
    !c     index is an integer array of dimension n
    !c       for i=1,...,nfree, index(i) are the indices of free variables
    !c       for i=nfree+1,...,n, index(i) are the indices of bound variables
    !c       on entry after the first iteration, index gives 
    !c         the free variables at the previous iteration.
    !c       on exit it gives the free variables based on the determination
    !c         in cauchy using the array iwhere.
    !c
    !c     indx2 is an integer array of dimension n
    !c       on entry indx2 is unspecified.
    !c       on exit with iter>0, indx2 indicates which variables
    !c          have changed status since the previous iteration.
    !c       for i= 1,...,nenter, indx2(i) have changed from bound to free.
    !c       for i= ileave+1,...,n, indx2(i) have changed from free to bound.
    !c 
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     ************
    integer :: iact,i,k

    nenter = 0
    ileave = n + 1
    if (iter > 0 .and. cnstnd) then
       !c                           count the entering and leaving variables.
       do i = 1, nfree
          k = indx(i)
          if (iwhere(k) > 0) then
             ileave = ileave - 1
             indx2(ileave) = k
             ! if (iprint >= 100) write (*,*) 'variable ',k,' leaves the set of free variables'
          end if
       end do
       do i = 1 + nfree, n
          k = indx(i)
          if (iwhere(k) <= 0) then
             nenter = nenter + 1
             indx2(nenter) = k
             ! if (iprint >= 100) write (*,*) 'variable ',k,' enters the set of free variables'
          end if
       end do
       ! if (iprint >= 99) write (*,*) n+1-ileave,' variables leave; ',nenter,' variables enter'
    end if
    wrk = (ileave < n+1) .or. (nenter > 0) .or. updatd

    !c     find the index set of free and active variables at the gcp.
    nfree = 0
    iact = n + 1    
    do i = 1, n
       if (iwhere(i) <= 0) then
          nfree = nfree + 1
          indx(nfree) = i
       else
          iact = iact - 1
          indx(iact) = i
       end if
    end do
    ! if (iprint >= 99) write (*,*) nfree,' variables are free at gcp ',iter + 1
    return
  end subroutine freev
  !c======================= The end of freev ==============================

  subroutine hpsolb(n, t, iorder, iheap)
    implicit none
    integer, intent(in) :: n
    real (dp), intent(inout) :: t(n)
    integer, intent(inout) :: iorder(n)
    integer, intent(in) :: iheap
    !c     ************
    !c
    !c     subroutine hpsolb 
    !c
    !c     this subroutine sorts out the least element of t, and puts the
    !c       remaining elements of t in a heap.
    !c 
    !c     n is an integer variable.
    !c       on entry n is the dimension of the arrays t and iorder.
    !c       on exit n is unchanged.
    !c
    !c     t is a double precision array of dimension n.
    !c       on entry t stores the elements to be sorted,
    !c       on exit t(n) stores the least elements of t, and t(1) to t(n-1)
    !c         stores the remaining elements in the form of a heap.
    !c
    !c     iorder is an integer array of dimension n.
    !c       on entry iorder(i) is the index of t(i).
    !c       on exit iorder(i) is still the index of t(i), but iorder may be
    !c         permuted in accordance with t.
    !c
    !c     iheap is an integer variable specifying the task.
    !c       on entry iheap should be set as follows:
    !c         iheap .eq. 0 if t(1) to t(n) is not in the form of a heap,
    !c         iheap .ne. 0 if otherwise.
    !c       on exit iheap is unchanged.
    !c
    !c
    !c     references:
    !c       algorithm 232 of cacm (j. w. j. williams): heapsort.
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c     ************
    integer :: i,j,k,indxin,indxou
    real (dp) :: ddum, output

    if (iheap == 0) then

       !c        rearrange the elements t(1) to t(n) to form a heap.

       do k = 2, n
          ddum = t(k)
          indxin = iorder(k)

          !c           add ddum to the heap.
          i = k
10        continue
          if (i > 1) then
             j = i/2
             if (ddum < t(j)) then
                t(i) = t(j)
                iorder(i) = iorder(j)
                i = j
                go to 10
             end if
          end if
          t(i) = ddum
          iorder(i) = indxin
       end do
    end if

    !c     assign to 'out' the value of t(1), the least member of the heap,
    !c        and rearrange the remaining members to form a heap as
    !c        elements 1 to n-1 of t.

    if (n > 1) then
       i = 1
       output = t(1)
       indxou = iorder(1)
       ddum = t(n)
       indxin = iorder(n)

       !c        restore the heap 
30     continue
       j = i+i
       if (j <= n-1) then
          if (t(j+1) < t(j)) j = j+1
          if (t(j) < ddum ) then
             t(i) = t(j)
             iorder(i) = iorder(j)
             i = j
             go to 30
          end if
       end if
       t(i) = ddum
       iorder(i) = indxin

       !c     put the least member in t(n).

       t(n) = output
       iorder(n) = indxou
    end if
    return
  end subroutine hpsolb
  !c====================== The end of hpsolb ==============================

  subroutine lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, &
       z, stp, dnorm, dtd,  xstep, stpmx, iter, ifun,&
       iback, nfgv, info, task, boxed, cnstnd, csave,&
       isave,dsave)
    implicit none
    integer, intent(in) :: n, nbd(n)
    real (dp), intent(in) :: l(n), u(n), f
    real (dp), intent(inout) :: x(n), g(n), d(n), r(n), t(n), z(n), dsave(13)
    real (dp), intent(inout) :: fold, gd, gdold, stp, dnorm, dtd, xstep, stpmx
    integer, intent(inout) :: iter, isave(2)
    integer, intent(inout) :: ifun, iback, nfgv, info
    character (len=60), intent(inout) :: task, csave
    logical, intent(inout) :: boxed, cnstnd
    !c     **********
    !c
    !c     subroutine lnsrlb
    !c
    !c     this subroutine calls subroutine dcsrch from the minpack2 library
    !c       to perform the line search.  subroutine dscrch is safeguarded so
    !c       that all trial points lie within the feasible region.
    !c
    !c     subprograms called:
    !c
    !c       minpack2 library ... dcsrch.
    !c
    !c       linpack ... dtrsl, ddot.
    !c
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     **********
    integer :: i
    real (dp) :: a1,a2
    real (dp), parameter :: one=1.0d0, zero=0.0d0, big=1.0d+10
    real (dp), parameter ::  ftol=1.0d-3, gtol=0.9d0, xtol=0.1d0


    if (task(1:5) == 'fg_ln') go to 556

    dtd = ddot(n,d,1,d,1)
    dnorm = sqrt(dtd)

    !c     determine the maximum step length.
    stpmx = big
    if (cnstnd) then
       if (iter == 0) then
          stpmx = one
       else
          do i = 1, n
             a1 = d(i)
             if (nbd(i) /= 0) then
                if (a1 < zero .and. nbd(i) <= 2) then
                   a2 = l(i) - x(i)
                   if (a2 >= zero) then
                      stpmx = zero
                   else if (a1*stpmx < a2) then
                      stpmx = a2/a1
                   end if
                else if (a1 > zero .and. nbd(i) >= 2) then
                   a2 = u(i) - x(i)
                   if (a2 <= zero) then
                      stpmx = zero
                   else if (a1*stpmx > a2) then
                      stpmx = a2/a1
                   end if
                end if
             end if
          end do
       end if
    end if

    if (iter == 0 .and. .not. boxed) then
       stp = min(one/dnorm, stpmx)
    else
       stp = one
    end if

    call dcopy(n,x,1,t,1)   
    call dcopy(n,g,1,r,1)    
    fold = f
    ifun = 0
    iback = 0
    csave = 'start'

556 continue
    gd = ddot(n,g,1,d,1)
    if (ifun == 0) then
       gdold=gd
       if (gd >= zero) then
          !c   the directional derivative >=0.
          !c   line search is impossible.
          info = -4
          return
       end if
    end if

    call dcsrch(f,gd,stp,ftol,gtol,xtol,zero,stpmx,csave,isave,dsave)    

    xstep = stp*dnorm

    if (csave(1:4) /= 'conv' .and. csave(1:4) /= 'warn') then
       task = 'fg_lnsrch'
       ifun = ifun + 1
       nfgv = nfgv + 1
       iback = ifun - 1
       if (stp == one) then
          call dcopy(n,z,1,x,1)
       else
          do i = 1, n
             x(i) = stp*d(i) + t(i)
          end do
       end if
    else
       task = 'new_x'
    end if

    return
  end subroutine lnsrlb
  !c======================= The end of lnsrlb =============================

  subroutine matupd(n, m, ws, wy, sy, ss, d, r, itail,&
       iupdat, col, head, theta, rr, dr, stp, dtd)
    implicit none
    integer, intent(in) :: n, m, iupdat
    real (dp), intent(inout) :: ws(n, m), wy(n, m), d(n), r(n)
    real (dp), intent(inout) :: sy(m, m), ss(m, m), theta
    integer, intent(inout) :: col, head, itail
    real (dp), intent(in) :: rr, dr, stp, dtd
    !c     ************
    !c
    !c     subroutine matupd
    !c
    !c       this subroutine updates matrices ws and wy, and forms the
    !c         middle matrix in b.
    !c
    !c     subprograms called:
    !c
    !c       linpack ... dcopy, ddot.
    !c
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     ************
    integer :: j,pointr
    real (dp), parameter :: one=1.0d0

    !c     set pointers for matrices ws and wy.

    if (iupdat <= m) then
       col = iupdat
       itail = mod(head+iupdat-2,m) + 1
    else
       itail = mod(itail,m) + 1
       head = mod(head,m) + 1
    end if

    !c     update matrices ws and wy.

    call dcopy(n,d,1,ws(1,itail),1)
    call dcopy(n,r,1,wy(1,itail),1)

    !c     set theta=yy/ys.

    theta = rr/dr

    !c     form the middle matrix in b.

    !c    update the upper triangle of ss,
    !c             and the lower triangle of sy:
    if (iupdat > m) then
       !c     move old information
       do j = 1, col - 1
          call dcopy(j,ss(2,j+1),1,ss(1,j),1)
          call dcopy(col-j,sy(j+1,j+1),1,sy(j,j),1)
       end do
    end if
    !c        add new information: the last row of sy
    !c                              and the last column of ss:
    pointr = head
    do j = 1, col - 1
       sy(col,j) = ddot(n,d,1,wy(1,pointr),1)
       ss(j,col) = ddot(n,ws(1,pointr),1,d,1)
       pointr = mod(pointr,m) + 1
    end do
    if (stp == one) then
       ss(col,col) = dtd
    else
       ss(col,col) = stp*stp*dtd
    end if
    sy(col,col) = dr
    return
  end subroutine matupd
  !c======================= The end of matupd =============================

  subroutine prn1lb(n, m, l, u, x, iprint, itfile, epsmch)
    implicit none
    integer, intent(in) :: n, m, iprint, itfile
    real(dp), intent(in) :: epsmch, x(n), l(n), u(n)
    !c     ************
    !c
    !c     subroutine prn1lb
    !c
    !c     this subroutine prints the input data, initial point, upper and
    !c       lower bounds of each variable, machine precision, as well as 
    !c       the headings of the output.
    !c
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     ************
    integer i

    !c     ************
    ! dummy, to suppress the message that the variable is not used
    i = itfile
    !c     ************

    if (iprint >= 0) then
       call dblepr("N, M and machine precision", -1,&
            (/dble(n),dble(m),epsmch/), 3)
       if(iprint >= 100) then
          call dblepr("L = ", -1, l, n)
          call dblepr("x0 = ", -1, x, n)
          call dblepr("u = ", -1, u, n)          
       end if
!!$       write (*,7001) epsmch
!!$       write (*,*) 'n = ',n,'    m = ',m
!!$       if (iprint >= 1) then
!!$          write (itfile,2001) epsmch
!!$          write (itfile,*)'n = ',n,'    m = ',m
!!$          write (itfile,9001)
!!$          if (iprint > 100) then
!!$             write (*,1004) 'l =',(l(i),i = 1,n)
!!$             write (*,1004) 'x0 =',(x(i),i = 1,n)
!!$             write (*,1004) 'u =',(u(i),i = 1,n)
!!$          endif
!!$       endif
    endif

!!$1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
!!$2001 format ('running the l-bfgs-b code',/,/,&
!!$         'it    = iteration number',/,&
!!$         'nf    = number of function evaluations',/,&
!!$         'nseg  = number of segments explored during the cauchy search',/,&
!!$         'nact  = number of active bounds at the generalized cauchy point'&
!!$         ,/,&
!!$         'sub   = manner in which the subspace minimization terminated:'&
!!$         ,/,'        con = converged, bnd = a bound was reached',/,&
!!$         'itls  = number of iterations performed in the line search',/,&
!!$         'stepl = step length used',/,&
!!$         'tstep = norm of the displacement (total step)',/,&
!!$         'projg = norm of the projected gradient',/,&
!!$         'f     = function value',/,/,&
!!$         '           * * *',/,/,&
!!$         'machine precision =',1p,d10.3)
!!$7001 format ('running the l-bfgs-b code',/,/,&
!!$          '           * * *',/,/,&
!!$          'machine precision =',1p,d10.3)
!!$9001 format (/,3x,'it',3x,'nf',2x,'nseg',2x,'nact',2x,'sub',2x,'itls',&
!!$          2x,'stepl',4x,'tstep',5x,'projg',8x,'f')

    return
  end subroutine prn1lb
  !c======================= The end of prn1lb =============================

  subroutine prn2lb(n, x, f, g, iprint, itfile, iter, nfgv, nact, &
       sbgnrm, nseg,   word, iword, iback, stp, xstep)
    implicit none
    integer, intent(in) :: n, iter, nfgv, nact
    real (dp), intent(in) :: x(n), f, g(n), sbgnrm, stp, xstep
    integer, intent(in) :: iprint, nseg, iword
    integer, intent(inout) :: itfile
    character (len=3), intent(out) :: word
    integer, intent(inout) :: iback
    !c     ************
    !c
    !c     subroutine prn2lb
    !c
    !c     this subroutine prints out new information after a successful
    !c       line search. 
    !c
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     ************
    integer :: imod, i
    real (dp) :: dm

    !c     ************
    ! dummy, to suppress the message that the variable is not used
    i = itfile
    i = nact
    i = nfgv
    i = nseg
    dm = stp    
    !c     ************

    !c           'word' records the status of subspace solutions.
    if (iword == 0) then
       !c     the subspace minimization converged.
       word = 'con'
    else if (iword == 1) then
       !c     the subspace minimization stopped at a bound.
       word = 'bnd'
    else if (iword == 5) then
       !c      the truncated newton step has been used.
       word = 'tnt'
    else
       word = '---'
    end if
    if (iprint >= 99) then
       call dblepr("line search times; norm of step", -1, &
            (/dble(iback), xstep/), 2)  
!!$       write (*,*) 'line search',iback,' times; norm of step = ',xstep
       call dblepr("at iterate /   f=    / |proj g|= ", &
            -1, (/dble(iter), f, sbgnrm/), 3)
       call dblepr("x = ", -1, x, n)
       call dblepr("g = ", -1, g, n)         
!!$       write (*,2001) iter,f,sbgnrm
!!$       if (iprint > 100) then
!!$          write (*,1004) 'x =',(x(i), i = 1, n)
!!$          write (*,1004) 'g =',(g(i), i = 1, n)
!!$       end if
    else if (iprint > 0) then
       imod = mod(iter,iprint)
!!$       if (imod == 0) write (*,2001) iter,f,sbgnrm
       if (imod == 0) call dblepr("at iterate /   f=    / |proj g|= ", &
            -1, (/dble(iter), f, sbgnrm/), 3)
    end if
!!$    if (iprint >= 1) write (itfile,3001)&
!!$         iter,nfgv,nseg,nact,word,iback,stp,xstep,sbgnrm,f

!!$1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
!!$2001 format ('at iterate',i5,4x,'f= ',f12.4,4x,'|proj g|= ',f12.4)
!!$3001 format(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),1p,2(1x,d10.3))
    return
  end subroutine prn2lb
  !c======================= The end of prn2lb =============================

  subroutine prn3lb(n, x, f, task, iprint, info, itfile,&
       iter, nfgv, nintol, nskip, nact, sbgnrm, &
       time, nseg, word, iback, stp, xstep, k, &
       cachyt, sbtime, lnscht)
    implicit none
    integer, intent(in) :: n, nact, nskip, nintol, nfgv,iter, iprint, info, nseg,iback, k
    real (dp), intent(in) :: x(n), f, time,sbgnrm, cachyt, stp, xstep, sbtime, lnscht
    character (len=60), intent(inout) :: task
    integer, intent(inout) :: itfile
    character (len=3), intent(in) :: word
    !c     ************
    !c
    !c     subroutine prn3lb
    !c
    !c     this subroutine prints out information when either a built-in
    !c       convergence test is satisfied or when an error message is
    !c       generated.
    !c       
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     ************

    if (task(1:5) == 'error') go to 999

    if (iprint >= 0) then
       call labelpr(' * * *', -1)
       call intpr1('tit = total number of iterations', -1, iter)
       call intpr1('tnf = total number of function evaluations', -1, nfgv)
       call intpr1('tnpt_quad = total number of segments explored during'// &
            ' cauchy searches', -1, nintol)
       call intpr1('skip = number of bfgs updates skipped', -1, nskip)
       call intpr1('nact = number of active bounds at final generalized'//&
            ' cauchy point', -1, nact)
       call dblepr1('projg = norm of the final projected gradient', -1, sbgnrm)
       call dblepr1('f = final function value', -1, f)
       call labelpr(' * * *', -1)
!!$       write (*,3003)
!!$       write (*,3004)
!!$       write (*,3005) n,iter,nfgv,nintol,nskip,nact,sbgnrm,f
       if (iprint >= 100) then
          call dblepr('x =', -1, x, n)
!!$          write (*,1004) 'x =',(x(i),i = 1,n)
       end if
!!$       if (iprint >= 1) write (*,*) ' f =',f
       if (iprint >= 1) call dblepr1('f =', -1, f)
    end if

999 continue

    if (iprint >= 0) then
       call labelpr(task, -1)
       select case (info)
       case(-1)
          call labelpr("Matrix in 1st Cholesky factorization in formk is not Pos. Def.",-1)
       case(-2)
          call labelpr("Matrix in 2st Cholesky factorization in formk is not Pos. Def.",-1)
       case(-3)
          call labelpr("Matrix in the Cholesky factorization in formt is not Pos. Def.",-1)
       case(-4)
          call labelpr("Derivative >= 0, backtracking line search impossible."//&
               ' previous x, f and g restored.'// &
               ' possible causes: 1 error in function or gradient evaluation;'// &
               ' 2 rounding errors dominate computation.',-1)
       case(-5)
          call labelpr("Warning:  more than 10 function and gradient evaluations"//&
               " in the last line search. termination"//&
               ' may possibly be caused by a bad search direction.', -1)
       case(-6)
          call intpr1("Input nbd(k) is invalid. k = ", -1, k)
       case(-7)
          call intpr1("l(k) > u(k).  No feasible solution. k =", -1, k)
       case(-8)
          call labelpr("The triangular system is singular.", -1)
       case(-9)
          call labelpr("Line search cannot locate an adequate point after"//&
               "  20 function and gradient evaluations"//&
               ' previous x, f and g restored.'// &
               ' possible causes: 1 error in function or gradient evaluation;'// &
               ' 2 rounding errors dominate computation.',-1)
       end select
!!$       write (*,3009) task
!!$       if (info /= 0) then
!!$          if (info == -1) write (*,9011)
!!$          if (info == -2) write (*,9012)
!!$          if (info == -3) write (*,9013)
!!$          if (info == -4) write (*,9014)
!!$          if (info == -5) write (*,9015)
!!$          if (info == -6) write (*,*)' input nbd(',k,') is invalid.'
!!$          if (info == -7) write (*,*)' l(',k,') > u(',k,').  no feasible solution.'
!!$          if (info == -8) write (*,9018)
!!$          if (info == -9) write (*,9019)
!!$       end if

       if (iprint >= 1) then
          !write (*,3007) cachyt,sbtime,lnscht
          call dblepr1(' cauchy time (seconds)', -1, cachyt)
          call dblepr1(' subspace minimization time (seconds)', -1, sbtime)
          call dblepr1(' line search time (seconds)', -1, lnscht)              
       end if

!!$       write (*,3008) time
       call dblepr1(' total user time (seconds)', -1, time)

!!$       if (iprint >= 1) then
!!$          if (info == -4 .or. info == -9) then
!!$             write (itfile,3002) iter,nfgv,nseg,nact,word,iback,stp,xstep
!!$          endif
!!$          write (itfile,3009) task
!!$          if (info /= 0) then
!!$             if (info .eq. -1) write (itfile,9011)
!!$             if (info .eq. -2) write (itfile,9012)
!!$             if (info .eq. -3) write (itfile,9013)
!!$             if (info .eq. -4) write (itfile,9014)
!!$             if (info .eq. -5) write (itfile,9015)
!!$             if (info .eq. -8) write (itfile,9018)
!!$             if (info .eq. -9) write (itfile,9019)
!!$          endif
!!$          write (itfile,3008) time
!!$       endif
    endif

!!$1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
!!$3002 format(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),6x,'-',10x,'-')
!!$3003 format (/, ' * * *',/,/, &
!!$         'tit = total number of iterations',/, &
!!$         'tnf = total number of function evaluations',/, &
!!$         'tnpt_quad = total number of segments explored during', ' cauchy searches',/, &
!!$         'skip = number of bfgs updates skipped',/, &
!!$         'nact = number of active bounds at final generalized', ' cauchy point',/, &
!!$         'projg = norm of the final projected gradient',/, &
!!$         'f = final function value',/,/, ' * * *')
!!$3004 format (/,3x,'n',3x,'tit',2x,'tnf',2x,'tnpt_quad',2x, &
!!$         'skip',2x,'nact',5x,'projg',8x,'f')
!!$3005 format (i5,2(1x,i4),(1x,i6),(2x,i4),(1x,i5),1p,2(2x,d10.3))
!!$3006 format (i5,2(1x,i4),2(1x,i6),(1x,i4),(1x,i5),7x,'-',10x,'-')
!!$3007 format (/,' cauchy time',1p,e10.3,' seconds.',/ &
!!$         ' subspace minimization time',1p,e10.3,' seconds.',/ &
!!$         ' line search time',1p,e10.3,' seconds.')
!!$3008 format (/,' total user time',1p,e10.3,' seconds.',/)
!!$3009 format (a60)
!!$9011 format (/, &
!!$         ' matrix in 1st cholesky factorization in formk is not pos. def.')
!!$9012 format (/, &
!!$         ' matrix in 2st cholesky factorization in formk is not pos. def.')
!!$9013 format (/, &
!!$         ' matrix in the cholesky factorization in formt is not pos. def.')
!!$9014 format (/, ' derivative >= 0, backtracking line search impossible.',/, &
!!$         ' previous x, f and g restored.',/, &
!!$         ' possible causes: 1 error in function or gradient evaluation;',/, &
!!$         ' 2 rounding errors dominate computation.')
!!$9015 format (/, ' warning: more than 10 function and gradient',/, &
!!$         ' evaluations in the last line search.  termination',/, &
!!$         ' may possibly be caused by a bad search direction.')
!!$9018 format (/,' the triangular system is singular.')
!!$9019 format (/, &
!!$         ' line search cannot locate an adequate point after 20 function',/ &
!!$         ,' and gradient evaluations.  previous x, f and g restored.',/, &
!!$         ' possible causes: 1 error in function or gradient evaluation;',/, &
!!$         ' 2 rounding error dominate computation.')
    return
  end subroutine prn3lb
  !c======================= The end of prn3lb =============================

  subroutine projgr(n, l, u, nbd, x, g, sbgnrm)
    implicit none
    integer, intent(in) :: n
    real (dp), intent(in) :: l(n), u(n), g(n)
    integer, intent(in) :: nbd(n)
    real (dp), intent(inout) :: x(n)
    real (dp), intent(inout) :: sbgnrm
    !c     ************
    !c
    !c     subroutine projgr
    !c
    !c     this subroutine computes the infinity norm of the projected
    !c       gradient.
    !c
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     ************
    integer :: i
    real (dp) :: gi
    real (dp), parameter :: one=1.0d0, zero=0.0d0

    sbgnrm = zero
    do i = 1, n
       gi = g(i)
       if (nbd(i) /= 0) then
          if (gi < zero) then
             if (nbd(i) >= 2) gi = max((x(i)-u(i)),gi)
          else
             if (nbd(i) <= 2) gi = min((x(i)-l(i)),gi)
          end if
       end if
       sbgnrm = max(sbgnrm,abs(gi))
    end do
    return
  end subroutine projgr
  !c======================= The end of projgr =============================

  subroutine subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, ws, wy,&
       theta, xx, gg, &
       col,  head, iword, wv,wn, iprint, info)
    implicit none
    integer, intent(in) :: n, m, nsub, ind(nsub), nbd(n), col, head, iprint
    real (dp), intent(in) :: l(n), u(n), ws(n, m), wy(n, m), theta
    real(dp), intent(in) :: xx(n), gg(n)
    real (dp), intent(inout) :: x(n), d(n), wn(2*m, 2*m), xp(n)
    integer, intent(inout) :: iword
    real (dp), intent(inout) :: wv(2*m)
    integer, intent(inout) ::  info
    !c ***************************************************************
    !c
    !c     this routine contains the major changes in the updated version.
    !c     the changes are described in the accompanying paper
    !c
    !c      jose luis morales, jorge nocedal
    !c      "remark on algorithm 788: l-bfgs-b: fortran subroutines for large-scale
    !c       bound constrained optimization". decemmber 27, 2010.
    !c
    !c             j.l. morales  departamento de matematicas, 
    !c                           instituto tecnologico autonomo de mexico
    !c                           mexico d.f.
    !c
    !c             j, nocedal    department of electrical engineering and
    !c                           computer science.
    !c                           northwestern university. evanston, il. usa
    !c
    !c                           january 17, 2011
    !c
    !c      ************************************************************
    !c                           
    !c
    !c     subroutine subsm
    !c
    !c     given xcp, l, u, r, an index set that specifies
    !c       the active set at xcp, and an l-bfgs matrix b 
    !c       (in terms of wy, ws, sy, wt, head, col, and theta), 
    !c       this subroutine computes an approximate solution
    !c       of the subspace problem
    !c
    !c       (p)   min q(x) = r'(x-xcp) + 1/2 (x-xcp)' b (x-xcp)
    !c
    !c             subject to l<=x<=u
    !c                       x_i=xcp_i for all i in a(xcp)
    !c                     
    !c       along the subspace unconstrained newton direction 
    !c       
    !c          d = -(z'bz)^(-1) r.
    !c
    !c       the formula for the newton direction, given the l-bfgs matrix
    !c       and the sherman-morrison formula, is
    !c
    !c          d = (1/theta)r + (1/theta*2) z'wk^(-1)w'z r.
    !c 
    !c       where
    !c                 k = [-d -y'zz'y/theta     l_a'-r_z'  ]
    !c                     [l_a -r_z           theta*s'aa's ]
    !c
    !c     note that this procedure for computing d differs 
    !c     from that described in [1]. one can show that the matrix k is
    !c     equal to the matrix m^[-1]n in that paper.
    !c
    !c     n is an integer variable.
    !c       on entry n is the dimension of the problem.
    !c       on exit n is unchanged.
    !c
    !c     m is an integer variable.
    !c       on entry m is the maximum number of variable metric corrections
    !c         used to define the limited memory matrix.
    !c       on exit m is unchanged.
    !c
    !c     nsub is an integer variable.
    !c       on entry nsub is the number of free variables.
    !c       on exit nsub is unchanged.
    !c
    !c     ind is an integer array of dimension nsub.
    !c       on entry ind specifies the coordinate indices of free variables.
    !c       on exit ind is unchanged.
    !c
    !c     l is a double precision array of dimension n.
    !c       on entry l is the lower bound of x.
    !c       on exit l is unchanged.
    !c
    !c     u is a double precision array of dimension n.
    !c       on entry u is the upper bound of x.
    !c       on exit u is unchanged.
    !c
    !c     nbd is a integer array of dimension n.
    !c       on entry nbd represents the type of bounds imposed on the
    !c         variables, and must be specified as follows:
    !c         nbd(i)=0 if x(i) is unbounded,
    !c                1 if x(i) has only a lower bound,
    !c                2 if x(i) has both lower and upper bounds, and
    !c                3 if x(i) has only an upper bound.
    !c       on exit nbd is unchanged.
    !c
    !c     x is a double precision array of dimension n.
    !c       on entry x specifies the cauchy point xcp. 
    !c       on exit x(i) is the minimizer of q over the subspace of
    !c                                                        free variables. 
    !c
    !c     d is a double precision array of dimension n.
    !c       on entry d is the reduced gradient of q at xcp.
    !c       on exit d is the newton direction of q. 
    !c
    !c    xp is a double precision array of dimension n.
    !c       used to safeguard the projected newton direction 
    !c
    !c    xx is a double precision array of dimension n
    !c       on entry it holds the current iterate
    !c       on output it is unchanged
    !c
    !c    gg is a double precision array of dimension n
    !c       on entry it holds the gradient at the current iterate
    !c       on output it is unchanged
    !c
    !c     ws and wy are double precision arrays;
    !c     theta is a double precision variable;
    !c     col is an integer variable;
    !c     head is an integer variable.
    !c       on entry they store the information defining the
    !c                                          limited memory bfgs matrix:
    !c         ws(n,m) stores s, a set of s-vectors;
    !c         wy(n,m) stores y, a set of y-vectors;
    !c         theta is the scaling factor specifying b_0 = theta i;
    !c         col is the number of variable metric corrections stored;
    !c         head is the location of the 1st s- (or y-) vector in s (or y).
    !c       on exit they are unchanged.
    !c
    !c     iword is an integer variable.
    !c       on entry iword is unspecified.
    !c       on exit iword specifies the status of the subspace solution.
    !c         iword = 0 if the solution is in the box,
    !c                 1 if some bound is encountered.
    !c
    !c     wv is a double precision working array of dimension 2m.
    !c
    !c     wn is a double precision array of dimension 2m x 2m.
    !c       on entry the upper triangle of wn stores the lel^t factorization
    !c         of the indefinite matrix
    !c
    !c              k = [-d -y'zz'y/theta     l_a'-r_z'  ]
    !c                  [l_a -r_z           theta*s'aa's ]
    !c                                                    where e = [-i  0]
    !c                                                              [ 0  i]
    !c       on exit wn is unchanged.
    !c
    !c     iprint is an integer variable that must be set by the user.
    !c       it controls the frequency and type of output generated:
    !c        iprint<0    no output is generated;
    !c        iprint=0    print only one line at the last iteration;
    !c        0<iprint<99 print also f and |proj g| every iprint iterations;
    !c        iprint=99   print details of every iteration except n-vectors;
    !c        iprint=100  print also the changes of active set and final x;
    !c        iprint>100  print details of every iteration including x and g;
    !c       when iprint > 0, the file iterate.dat will be created to
    !c                        summarize the iteration.
    !c
    !c     info is an integer variable.
    !c       on entry info is unspecified.
    !c       on exit info = 0       for normal return,
    !c                    = nonzero for abnormal return 
    !c                                  when the matrix k is ill-conditioned.
    !c
    !c     subprograms called:
    !c
    !c       linpack dtrsl.
    !c
    !c
    !c     references:
    !c
    !c       [1] r. h. byrd, p. lu, j. nocedal and c. zhu, ``a limited
    !c       memory algorithm for bound constrained optimization'',
    !c       siam j. scientifi!c     omputing 16 (1995), no. 5, pp. 1190--1208.
    !c
    !c
    !c
    !c                           *  *  *
    !c
    !c     neos, november 1994. (latest revision june 1996.)
    !c     optimization technology center.
    !c     argonne national laboratory and northwestern university.
    !c     written by
    !c                        ciyou zhu
    !c     in collaboration with r.h. byrd, p. lu-chen and j. nocedal.
    !c
    !c
    !c     ************
    integer :: pointr,m2,col2,ibd,jy,js,i,j,k
    real (dp) :: alpha,xk,dk,temp1,temp2, dd_p
    real (dp), parameter :: one=1.0d0, zero=0.0d0

    if (nsub <= 0) return
    ! if (iprint >= 99) write (*,1001)
    if (iprint >= 99) call labelpr("----------------subsm entered-----------------",-1)

    !c     compute wv = w'zd.

    pointr = head
    do i = 1, col
       temp1 = zero
       temp2 = zero
       do j = 1, nsub
          k = ind(j)
          temp1 = temp1 + wy(k,pointr)*d(j)
          temp2 = temp2 + ws(k,pointr)*d(j)
       end do
       wv(i) = temp1
       wv(col + i) = theta*temp2
       pointr = mod(pointr,m) + 1
    end do

    !c     compute wv:=k^(-1)wv.

    m2 = 2*m
    col2 = 2*col
    call dtrsl(wn,m2,col2,wv,11,info)
    if (info /= 0) return
    do i = 1, col
       wv(i) = -wv(i)
    end do
    call dtrsl(wn,m2,col2,wv,01,info)
    if (info /= 0) return

    !c     compute d = (1/theta)d + (1/theta**2)z'w wv.

    pointr = head
    do jy = 1, col
       js = col + jy
       do i = 1, nsub
          k = ind(i)
          d(i) = d(i) + wy(k,pointr)*wv(jy)/theta + ws(k,pointr)*wv(js)
       end do
       pointr = mod(pointr,m) + 1
    end do

    ! d = d*1/theta
    call dscal( nsub, one/theta, d, 1 )

    goto 12345
    !------------------------------------------------------------------
    !
    !                      Added in version 3.0  (BEGIN)
    !
    !  R 4.0 does not have this part.
    !  I am ingoring this part because the code does not behave well
    !  when I include it.  (april, 2020)
    !
    !------------------------------------------------------------------
    !c 
    !c-----------------------------------------------------------------
    !c     let us try the projection, d is the newton direction
    iword = 0

    call dcopy ( n, x, 1, xp, 1 )


    do i=1, nsub
       k  = ind(i)
       dk = d(i)
       xk = x(k)
       if ( nbd(k) /= 0 ) then
          !c
          if ( nbd(k) == 1 ) then          ! lower bounds only
             x(k) = max( l(k), xk + dk )
             if ( x(k) == l(k) ) iword = 1
          else 
             !c     
             if ( nbd(k) == 2 ) then       ! upper and lower bounds
                xk   = max( l(k), xk + dk ) 
                x(k) = min( u(k), xk )
                if ( x(k) == l(k) .or. x(k) == u(k) ) iword = 1
             else
                !c
                if ( nbd(k) == 3 ) then    ! upper bounds only
                   x(k) = min( u(k), xk + dk )
                   if ( x(k) == u(k) ) iword = 1
                end if
             end if
          end if
          !c            
       else                                ! free variables
          x(k) = xk + dk
       end if
    end do
    !c
    if ( iword == 0 )   go to 911    
    !c
    !c     check sign of the directional derivative
    !c
    dd_p = zero
    do  i=1, n
       dd_p  = dd_p + (x(i) - xx(i))*gg(i)
    end do
    if ( dd_p > zero ) then
       call dcopy( n, xp, 1, x, 1 )
       ! write (*,*) ' positive dir derivative in projection '
       ! write (*,*) ' using the backtracking step '
    else
       go to 911
    endif

    !------------------------------------------------------------------
    !
    !                      Added in version 3.0 (END)
    !
    !------------------------------------------------------------------

12345 continue       

    !c
    !c-----------------------------------------------------------------
    !c
    alpha = one
    temp1 = alpha
    ibd = 0
    do i = 1, nsub
       k = ind(i)
       dk = d(i)
       if (nbd(k) /= 0) then
          if (dk < zero .and. nbd(k) <= 2) then
             temp2 = l(k) - x(k)
             if (temp2 >= zero) then
                temp1 = zero
             else if (dk*alpha < temp2) then
                temp1 = temp2/dk
             end if
          else if (dk > zero .and. nbd(k) >= 2) then
             temp2 = u(k) - x(k)
             if (temp2 <= zero) then
                temp1 = zero
             else if (dk*alpha > temp2) then
                temp1 = temp2/dk
             end if
          end if
          if (temp1 < alpha) then
             alpha = temp1
             ibd = i
          end if
       end if
    end do

    if (alpha < one) then
       dk = d(ibd)
       k = ind(ibd)
       if (dk > zero) then
          x(k) = u(k)
          d(ibd) = zero
       else if (dk < zero) then
          x(k) = l(k)
          d(ibd) = zero
       end if
    end if
    do i = 1, nsub
       k = ind(i)
       x(k) = x(k) + alpha*d(i)
    end do

    iword = 0
    if(alpha < one) iword = 1

911 continue

    ! if (iprint >= 99) write (*,1004)
    if (iprint >= 99) call labelpr("----------------exit subsm-----------------",-1)

    !1001 format (/,'----------------subsm entered-----------------',/)
    !1004 format (/,'----------------exit subsm --------------------',/)
    return
  end subroutine subsm
  !c====================== The end of subsm ===============================

  subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,&
       task,isave,dsave)    
    implicit none
    real (dp), intent(in) :: f, g, ftol, gtol, xtol, stpmin, stpmax
    real (dp), intent(inout) :: stp
    character (len=*), intent(inout) :: task
    integer, intent(inout) :: isave(2)
    real (dp), intent(inout) :: dsave(13)
    !c    **********
    !c
    !c    subroutine dcsrch
    !c
    !c    this subroutine finds a step that satisfies a sufficient
    !c    decrease condition and a curvature condition.
    !c
    !c    each call of the subroutine updates an interval with 
    !c    endpoints stx and sty. the interval is initially chosen 
    !c    so that it contains a minimizer of the modified function
    !c
    !c          psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
    !c
    !c    if psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
    !c    interval is chosen so that it contains a minimizer of f. 
    !c
    !c    the algorithm is designed to find a step that satisfies 
    !c    the sufficient decrease condition 
    !c
    !c          f(stp) <= f(0) + ftol*stp*f'(0),
    !c
    !c    and the curvature condition
    !c
    !c          abs(f'(stp)) <= gtol*abs(f'(0)).
    !c
    !c    if ftol is less than gtol and if, for example, the function
    !c    is bounded below, then there is always a step which satisfies
    !c    both conditions. 
    !c
    !c    if no step can be found that satisfies both conditions, then 
    !c    the algorithm stops with a warning. in this case stp only 
    !c    satisfies the sufficient decrease condition.
    !c
    !c    a typical invocation of dcsrch has the following outline:
    !c
    !c    task = 'start'
    !c  10 continue
    !c       call dcsrch( ... )
    !c       if (task .eq. 'fg') then
    !c          evaluate the function and the gradient at stp 
    !c          goto 10
    !c          end if
    !c
    !c    note: the user must no alter work arrays between calls.
    !c
    !c    the subroutine statement is
    !c
    !c       subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
    !c                         task,isave,dsave)
    !c    where
    !c
    !c      f is a double precision variable.
    !c        on initial entry f is the value of the function at 0.
    !c           on subsequent entries f is the value of the 
    !c           function at stp.
    !c        on exit f is the value of the function at stp.
    !c
    !c      g is a double precision variable.
    !c        on initial entry g is the derivative of the function at 0.
    !c           on subsequent entries g is the derivative of the 
    !c           function at stp.
    !c        on exit g is the derivative of the function at stp.
    !c
    !c      stp is a double precision variable. 
    !c        on entry stp is the current estimate of a satisfactory 
    !c           step. on initial entry, a positive initial estimate 
    !c           must be provided. 
    !c        on exit stp is the current estimate of a satisfactory step
    !c           if task = 'fg'. if task = 'conv' then stp satisfies
    !c           the sufficient decrease and curvature condition.
    !c
    !c      ftol is a double precision variable.
    !c        on entry ftol specifies a nonnegative tolerance for the 
    !c           sufficient decrease condition.
    !c        on exit ftol is unchanged.
    !c
    !c      gtol is a double precision variable.
    !c        on entry gtol specifies a nonnegative tolerance for the 
    !c           curvature condition. 
    !c        on exit gtol is unchanged.
    !c
    !c      xtol is a double precision variable.
    !c        on entry xtol specifies a nonnegative relative tolerance
    !c           for an acceptable step. the subroutine exits with a
    !c           warning if the relative difference between sty and stx
    !c           is less than xtol.
    !c        on exit xtol is unchanged.
    !c
    !c      stpmin is a double precision variable.
    !c        on entry stpmin is a nonnegative lower bound for the step.
    !c        on exit stpmin is unchanged.
    !c
    !c      stpmax is a double precision variable.
    !c        on entry stpmax is a nonnegative upper bound for the step.
    !c        on exit stpmax is unchanged.
    !c
    !c      task is a character variable of length at least 60.
    !c        on initial entry task must be set to 'start'.
    !c        on exit task indicates the required action:
    !c
    !c           if task(1:2) = 'fg' then evaluate the function and 
    !c           derivative at stp and call dcsrch again.
    !c
    !c           if task(1:4) = 'conv' then the search is successful.
    !c
    !c           if task(1:4) = 'warn' then the subroutine is not able
    !c           to satisfy the convergence conditions. the exit value of
    !c           stp contains the best point found during the search.
    !c
    !c           if task(1:5) = 'error' then there is an error in the
    !c           input arguments.
    !c
    !c        on exit with convergence, a warning or an error, the
    !c           variable task contains additional information.
    !c
    !c      isave is an integer work array of dimension 2.
    !c        
    !c      dsave is a double precision work array of dimension 13.
    !c
    !c    subprograms called
    !c
    !c      minpack-2 ... dcstep
    !c
    !c    minpack-1 project. june 1983.
    !c    argonne national laboratory. 
    !c    jorge j. more' and david j. thuente.
    !c
    !c    minpack-2 project. october 1993.
    !c    argonne national laboratory and university of minnesota. 
    !c    brett m. averick, richard g. carter, and jorge j. more'. 
    !c
    !c    **********
    real (dp), parameter :: zero=0.0d0, p5=0.5d0, p66=0.66d0, xtrapl=1.1d0, xtrapu=4.0d0
    logical :: brackt
    integer :: stage
    real (dp) :: finit,ftest,fm,fx,fxm,fy,fym,ginit,gtest
    real (dp) :: gm,gx,gxm,gy,gym,stx,sty,stmin,stmax,width,width1


    ! !c     initialization block.

    if (task(1:5) == 'start') then

       !c        check the input arguments for errors.

       if (stp < stpmin) task = 'error: stp .lt. stpmin'
       if (stp > stpmax) task = 'error: stp .gt. stpmax'
       if (g >= zero) task = 'error: initial g .ge. zero'
       if (ftol < zero) task = 'error: ftol .lt. zero'
       if (gtol < zero) task = 'error: gtol .lt. zero'
       if (xtol < zero) task = 'error: xtol .lt. zero'
       if (stpmin < zero) task = 'error: stpmin .lt. zero'
       if (stpmax < stpmin) task = 'error: stpmax .lt. stpmin'

       !c        exit if there are errors on input.

       if (task(1:5) == 'error') return

       !c        initialize local variables.

       brackt = .false.
       stage = 1
       finit = f
       ginit = g
       gtest = ftol*ginit
       width = stpmax - stpmin
       width1 = width/p5

       !c        the variables stx, fx, gx contain the values of the step, 
       !c        function, and derivative at the best step. 
       !c        the variables sty, fy, gy contain the value of the step, 
       !c        function, and derivative at sty.
       !c        the variables stp, f, g contain the values of the step, 
       !c        function, and derivative at stp.

       stx = zero
       fx = finit
       gx = ginit
       sty = zero
       fy = finit
       gy = ginit
       stmin = zero
       stmax = stp + xtrapu*stp
       task = 'fg'

       go to 1000

    else

       !c        restore local variables.

       if (isave(1) == 1) then
          brackt = .true.
       else
          brackt = .false.
       end if
       stage = isave(2)
       ginit = dsave(1)
       gtest = dsave(2)
       gx = dsave(3)
       gy = dsave(4)
       finit = dsave(5)
       fx = dsave(6)
       fy = dsave(7)
       stx = dsave(8)
       sty = dsave(9)
       stmin = dsave(10)
       stmax = dsave(11)
       width = dsave(12)
       width1 = dsave(13)

    end if

    !c     if psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
    !c     algorithm enters the second stage.

    ftest = finit + stp*gtest
    if (stage == 1 .and. f <= ftest .and. g >= zero) stage = 2


    !c     test for warnings.

    if (brackt .and. (stp <= stmin .or. stp >= stmax)) &
         task = 'warning: rounding errors prevent progress'
    if (brackt .and. stmax - stmin <= xtol*stmax) &
         task = 'warning: xtol test satisfied'
    if (stp == stpmax .and. f <= ftest .and. g <= gtest)&
         task = 'warning: stp = stpmax'
    if (stp == stpmin .and. (f > ftest .or. g >= gtest)) &
         task = 'warning: stp = stpmin'

    !c     test for convergence.
    if (f <= ftest .and. abs(g) <= gtol*(-ginit)) task = 'convergence'

    !c     test for termination.

    if (task(1:4) == 'warn' .or. task(1:4) == 'conv') go to 1000

    !c     a modified function is used to predict the step during the
    !c     first stage if a lower function value has been obtained but 
    !c     the decrease is not sufficient.

    if (stage == 1 .and. f <= fx .and. f > ftest) then

       !c        define the modified function and derivative values.

       fm = f - stp*gtest
       fxm = fx - stx*gtest
       fym = fy - sty*gtest
       gm = g - gtest
       gxm = gx - gtest
       gym = gy - gtest

       !c        call dcstep to update stx, sty, and to compute the new step.

       call dcstep(stx,fxm,gxm,sty,fym,gym,stp,fm,gm, brackt,stmin,stmax)

       !c        reset the function and derivative values for f.

       fx = fxm + stx*gtest
       fy = fym + sty*gtest
       gx = gxm + gtest
       gy = gym + gtest

    else

       !c       call dcstep to update stx, sty, and to compute the new step.

       call dcstep(stx,fx,gx,sty,fy,gy,stp,f,g, brackt,stmin,stmax)

    end if

    !c     decide if a bisection step is needed.

    if (brackt) then
       if (abs(sty-stx) >= p66*width1) stp = stx + p5*(sty - stx)
       width1 = width
       width = abs(sty-stx)
    end if

    !c     set the minimum and maximum steps allowed for stp.

    if (brackt) then
       stmin = min(stx,sty)
       stmax = max(stx,sty)
    else
       stmin = stp + xtrapl*(stp - stx)
       stmax = stp + xtrapu*(stp - stx)
    end if

    !c     force the step to be within the bounds stpmax and stpmin.

    stp = max(stp,stpmin)
    stp = min(stp,stpmax)

    !c     if further progress is not possible, let stp be the best
    !c     point obtained during the search.

    if (brackt .and. (stp <= stmin .or. stp >= stmax) &
         .or. (brackt .and. stmax-stmin <= xtol*stmax)) stp = stx

    !c     obtain another function and derivative.

    task = 'fg'

1000 continue

    !c     save local variables.

    if (brackt) then
       isave(1) = 1
    else
       isave(1) = 0
    end if
    isave(2) = stage
    dsave(1) = ginit
    dsave(2) = gtest
    dsave(3) = gx
    dsave(4) = gy
    dsave(5) = finit
    dsave(6) = fx
    dsave(7) = fy
    dsave(8) = stx
    dsave(9) = sty
    dsave(10) = stmin
    dsave(11) = stmax
    dsave(12) = width
    dsave(13) = width1
    return
  end subroutine dcsrch
  !c====================== The end of dcsrch ==============================


  subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dpvalue,brackt, stpmin,stpmax)
    implicit none
    real (dp), intent(inout) :: stx, fx, dx, dy, stp
    real (dp), intent(inout) :: sty, fy
    real (dp), intent(in) :: fp, dpvalue, stpmin, stpmax
    logical, intent(inout) :: brackt
    !c     **********
    !c
    !c     subroutine dcstep
    !c
    !c     this subroutine computes a safeguarded step for a search
    !c     procedure and updates an interval that contains a step that
    !c     satisfies a sufficient decrease and a curvature condition.
    !c
    !c     the parameter stx contains the step with the least function
    !c     value. if brackt is set to .true. then a minimizer has
    !c     been bracketed in an interval with endpoints stx and sty.
    !c     the parameter stp contains the current step. 
    !c     the subroutine assumes that if brackt is set to .true. then
    !c
    !c           min(stx,sty) < stp < max(stx,sty),
    !c
    !c     and that the derivative at stx is negative in the direction 
    !c     of the step.
    !c
    !c     the subroutine statement is
    !c
    !c       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
    !c                         stpmin,stpmax)
    !c
    !c     where
    !c
    !c       stx is a double precision variable.
    !c         on entry stx is the best step obtained so far and is an
    !c            endpoint of the interval that contains the minimizer. 
    !c         on exit stx is the updated best step.
    !c
    !c       fx is a double precision variable.
    !c         on entry fx is the function at stx.
    !c         on exit fx is the function at stx.
    !c
    !c       dx is a double precision variable.
    !c         on entry dx is the derivative of the function at 
    !c            stx. the derivative must be negative in the direction of 
    !c            the step, that is, dx and stp - stx must have opposite 
    !c            signs.
    !c         on exit dx is the derivative of the function at stx.
    !c
    !c       sty is a double precision variable.
    !c         on entry sty is the second endpoint of the interval that 
    !c            contains the minimizer.
    !c         on exit sty is the updated endpoint of the interval that 
    !c            contains the minimizer.
    !c
    !c       fy is a double precision variable.
    !c         on entry fy is the function at sty.
    !c         on exit fy is the function at sty.
    !c
    !c       dy is a double precision variable.
    !c         on entry dy is the derivative of the function at sty.
    !c         on exit dy is the derivative of the function at the exit sty.
    !c
    !c       stp is a double precision variable.
    !c         on entry stp is the current step. if brackt is set to .true.
    !c            then on input stp must be between stx and sty. 
    !c         on exit stp is a new trial step.
    !c
    !c       fp is a double precision variable.
    !c         on entry fp is the function at stp
    !c         on exit fp is unchanged.
    !c
    !c       dp is a double precision variable.
    !c         on entry dp is the the derivative of the function at stp.
    !c         on exit dp is unchanged.
    !c
    !c       brackt is an logical variable.
    !c         on entry brackt specifies if a minimizer has been bracketed.
    !c            initially brackt must be set to .false.
    !c         on exit brackt specifies if a minimizer has been bracketed.
    !c            when a minimizer is bracketed brackt is set to .true.
    !c
    !c       stpmin is a double precision variable.
    !c         on entry stpmin is a lower bound for the step.
    !c         on exit stpmin is unchanged.
    !c
    !c       stpmax is a double precision variable.
    !c         on entry stpmax is an upper bound for the step.
    !c         on exit stpmax is unchanged.
    !c
    !c     minpack-1 project. june 1983
    !c     argonne national laboratory. 
    !c     jorge j. more' and david j. thuente.
    !c
    !c     minpack-2 project. october 1993.
    !c     argonne national laboratory and university of minnesota. 
    !c     brett m. averick and jorge j. more'.
    !c
    !c     **********
    real (dp), parameter :: zero=0.0d0, p66=0.66d0,two=2.0d0,three=3.0d0
    real (dp) :: gama,p,q,r,s,sgnd,stpc,stpf,stpq,theta

    sgnd = dpvalue*(dx/abs(dx))

    !c     first case: a higher function value. the minimum is bracketed. 
    !c     if the cubic step is closer to stx than the quadratic step, the 
    !c     cubic step is taken, otherwise the average of the cubic and 
    !c     quadratic steps is taken.

    if (fp > fx) then
       theta = three*(fx - fp)/(stp - stx) + dx + dpvalue
       s = max(abs(theta),abs(dx),abs(dpvalue))
       gama = s*sqrt((theta/s)**2 - (dx/s)*(dpvalue/s))
       if (stp < stx) gama = -gama
       p = (gama - dx) + theta
       q = ((gama - dx) + gama) + dpvalue
       r = p/q
       stpc =stx + r*(stp - stx)
       stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/two)* (stp - stx)
       if (abs(stpc-stx) < abs(stpq-stx)) then
          stpf = stpc
       else
          stpf = stpc + (stpq - stpc)/two
       end if
       brackt = .true.


       !c     second case: a lower function value and derivatives of opposite 
       !c     sign. the minimum is bracketed. if the cubic step is farther from 
       !c     stp than the secant step, the cubic step is taken, otherwise the 
       !c     secant step is taken.

    else if (sgnd < zero) then
       theta = three*(fx - fp)/(stp - stx) + dx + dpvalue
       s = max(abs(theta),abs(dx),abs(dpvalue))
       gama = s*sqrt((theta/s)**2 - (dx/s)*(dpvalue/s))
       if (stp > stx) gama = -gama
       p = (gama - dpvalue) + theta
       q = ((gama - dpvalue) + gama) + dx
       r = p/q
       stpc= stp + r*(stx - stp)
       stpq = stp + (dpvalue/(dpvalue - dx))*(stx - stp)
       if (abs(stpc-stp) > abs(stpq-stp)) then
          stpf = stpc
       else
          stpf = stpq
       end if
       brackt = .true.

       !c     third case: a lower function value, derivatives of the same sign,
       !c     and the magnitude of the derivative decreases.

    else if (abs(dpvalue) < abs(dx)) then

       !c        the cubic step is computed only if the cubic tends to infinity 
       !c        in the direction of the step or if the minimum of the cubic
       !c        is beyond stp. otherwise the cubic step is defined to be the 
       !c        secant step.

       theta = three*(fx - fp)/(stp - stx) + dx + dpvalue
       s = max(abs(theta),abs(dx),abs(dpvalue))
       gama = s*sqrt(max(zero,(theta/s)**2-(dx/s)*(dpvalue/s)))

       !c        the case gamma = 0 only arises if the cubic does not tend
       !c        to infinity in the direction of the step.

       if (stp > stx) gama = -gama
       p = (gama - dpvalue) + theta
       q = (gama + (dx - dpvalue)) + gama
       r = p/q
       if (r < zero .and. gama /= zero) then
          stpc = stp + r*(stx - stp)
       else if (stp > stx) then
          stpc = stpmax
       else
          stpc = stpmin
       end if
       stpq = stp + (dpvalue/(dpvalue - dx))*(stx - stp)

       if (brackt) then

          !c           a minimizer has been bracketed. if the cubic step is 
          !c           closer to stp than the secant step, the cubic step is 
          !c           taken, otherwise the secant step is taken.

          if (abs(stpc-stp) < abs(stpq-stp)) then
             stpf = stpc
          else
             stpf = stpq
          end if
          if (stp > stx) then
             stpf = min(stp+p66*(sty-stp),stpf)
          else
             stpf = max(stp+p66*(sty-stp),stpf)
          end if
       else

          !c           a minimizer has not been bracketed. if the cubic step is 
          !c           farther from stp than the secant step, the cubic step is 
          !c           taken, otherwise the secant step is taken.

          if (abs(stpc-stp) > abs(stpq-stp)) then
             stpf = stpc
          else
             stpf = stpq
          end if
          stpf = min(stpmax,stpf)
          stpf = max(stpmin,stpf)
       end if

       !c     fourth case: a lower function value, derivatives of the same sign, 
       !c     and the magnitude of the derivative does not decrease. if the 
       !c     minimum is not bracketed, the step is either stpmin or stpmax, 
       !c     otherwise the cubic step is taken.

    else
       if (brackt) then
          theta = three*(fp - fy)/(sty - stp) + dy + dpvalue
          s = max(abs(theta),abs(dy),abs(dpvalue))
          gama = s*sqrt((theta/s)**2 - (dy/s)*(dpvalue/s))
          if (stp > sty) gama = -gama
          p = (gama - dpvalue) + theta
          q = ((gama - dpvalue) + gama) + dy
          r = p/q
          stpc = stp + r*(sty - stp)
          stpf = stpc
       else if (stp > stx) then
          stpf = stpmax
       else
          stpf = stpmin
       end if
    end if

    !c     update the interval which contains a minimizer.

    if (fp > fx) then
       sty = stp
       fy = fp
       dy = dpvalue
    else
       if (sgnd < zero) then
          sty = stx
          fy = fx
          dy = dx
       end if
       stx = stp
       fx = fp
       dx = dpvalue
    end if

    !c     compute the new step.

    stp = stpf
    return
  end subroutine dcstep
  !c====================== The end of dcstep ==============================

  subroutine timer(ttime)
    implicit none
    real (dp), intent(inout) :: ttime
    !c     this routine computes cpu time in double precision; it makes use of 
    !c     the intrinsic f90 cpu_time therefore a conversion type is
    !c     needed.
    !c
    !c           j.l morales  departamento de matematicas, 
    !c             instituto tecnologico autonomo de mexico
    !c             mexico d.f.
    !c
    !c           j.l nocedal  department of electrical engineering and
    !c             computer science.
    !c             northwestern university. evanston, il. usa
    !c              
    !c             january 21, 2011
    !c
    real :: temp      
    temp = sngl(ttime)
    call cpu_time(temp)
    ttime = dble(temp)
    return
  end subroutine timer


end module lbfgsb
