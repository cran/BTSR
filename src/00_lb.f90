module lpk

  implicit none
  public
  private :: dp

  integer, parameter :: dp = kind(1.d0)

contains 	
  !----------------------------------------------------------------
  !
  !   subroutines from linpack library 
  !
  !----------------------------------------------------------------
  subroutine dtrsl(t,ldt,n,b,job,info)
    ! c
    ! c
    ! c     dtrsl solves systems of the form
    ! c
    ! c                   t * x = b
    ! c     or
    ! c                   trans(t) * x = b
    ! c
    ! c     where t is a triangular matrix of order n. here trans(t)
    ! c     denotes the transpose of the matrix t.
    ! c
    ! c     on entry
    ! c
    ! c         t         double precision(ldt,n)
    ! c                   t contains the matrix of the system. the zero
    ! c                   elements of the matrix are not referenced, and
    ! c                   the corresponding elements of the array can be
    ! c                   used to store other information.
    ! c
    ! c         ldt       integer
    ! c                   ldt is the leading dimension of the array t.
    ! c
    ! c         n         integer
    ! c                   n is the order of the system.
    ! c
    ! c         b         double precision(n).
    ! c                   b contains the right hand side of the system.
    ! c
    ! c         job       integer
    ! c                   job specifies what kind of system is to be solved.
    ! c                   if job is
    ! c
    ! c                        00   solve t*x=b, t lower triangular,
    ! c                        01   solve t*x=b, t upper triangular,
    ! c                        10   solve trans(t)*x=b, t lower triangular,
    ! c                        11   solve trans(t)*x=b, t upper triangular.
    ! c
    ! c     on return
    ! c
    ! c         b         b contains the solution, if info .eq. 0.
    ! c                   otherwise b is unaltered.
    ! c
    ! c         info      integer
    ! c                   info contains zero if the system is nonsingular.
    ! c                   otherwise info contains the index of
    ! c                   the first zero diagonal element of t.
    ! c
    ! c     linpack. this version dated 08/14/78 .
    ! c     g. w. stewart, university of maryland, argonne national lab.
    ! c
    !----------------------------------------------------------------------------
    ! porto alegre. december 8 , 2012.
    ! changed b(1) for b(n)
    ! porto alegre. july 15, 2016.
    ! changed t(ldt,1) for t(ldt,n)
    !
    implicit none
    integer, intent(in) :: ldt, n, job
    real (dp), intent(in) :: t(ldt,n)
    real (dp), intent(inout) :: b(n)
    integer, intent(inout) :: info
    real (dp) :: temp
    integer :: cases,j,jj

    ! c
    ! c     begin block permitting ...exits to 150
    ! c
    ! c        check for zero diagonal elements.
    ! c
    do info = 1, n
       if (t(info,info) == 0.0d0) go to 150
    end do
    info = 0
    ! c
    ! c        determine the task and go to it.
    ! c      
    cases = 1
    if (mod(job,10) /= 0) cases = 2
    if (mod(job,100)/10 /= 0) cases = cases + 2
    select case ( cases )
    case ( 1)
       go to 20
    case ( 2)
       go to 50
    case ( 3)
       go to 80
    case ( 4)
       go to 110
    end select
    ! c
    ! c        solve t*x=b for t lower triangular
    ! c
20  continue
    b(1) = b(1)/t(1,1)
    if (n < 2) go to 40
    do j = 2, n
       temp = -b(j-1)
       call daxpy(n-j+1,temp,t(j,j-1),1,b(j),1)
       b(j) = b(j)/t(j,j)
    end do
40  continue
    go to 140
    ! c
    ! c        solve t*x=b for t upper triangular.
    ! c      
50  continue
    b(n) = b(n)/t(n,n)
    if (n < 2) go to 70
    do jj = 2, n
       j = n - jj + 1
       temp = -b(j+1)
       call daxpy(j,temp,t(1,j+1),1,b(1),1)
       b(j) = b(j)/t(j,j)
    end do
70  continue
    go to 140
    ! c
    ! c        solve trans(t)*x=b for t lower triangular.
    ! c      
80  continue
    b(n) = b(n)/t(n,n)
    if (n < 2) go to 100
    do jj = 2, n
       j = n - jj + 1
       b(j) = b(j) - ddot(jj-1,t(j+1,j),1,b(j+1),1)
       b(j) = b(j)/t(j,j)
    end do
100 continue
    go to 140
    ! c
    ! c        solve trans(t)*x=b for t upper triangular.
    ! c      
110 continue
    b(1) = b(1)/t(1,1)
    if (n < 2) go to 130
    do j = 2, n
       b(j) = b(j) - ddot(j-1,t(1,j),1,b(1),1)
       b(j) = b(j)/t(j,j)
    end do
130 continue
140 continue
150 continue
    return
  end subroutine dtrsl

  subroutine dpofa(a,lda,n, info)
    ! Porto Alegre. December 8 , 2012.
    ! changed a(lda,1) for a(lda,n).
    ! Taiane Schaedler Prass
    !-------------------------------------------------------------------
    ! c
    ! c     dpofa factors a double precision symmetric positive definite
    ! c     matrix.
    ! c
    ! c     dpofa is usually called by dpoco, but it can be called
    ! c     directly with a saving in time if  rcond  is not needed.
    ! c     (time for dpoco) = (1 + 18/n)*(time for dpofa) .
    ! c
    ! c     on entry
    ! c
    ! c        a       double precision(lda, n)
    ! c                the symmetric matrix to be factored.  only the
    ! c                diagonal and upper triangle are used.
    ! c
    ! c        lda     integer
    ! c                the leading dimension of the array  a .
    ! c
    ! c        n       integer
    ! c                the order of the matrix  a .
    ! c
    ! c     on return
    ! c
    ! c        a       an upper triangular matrix  r  so that  a = trans(r)*r
    ! c                where  trans(r)  is the transpose.
    ! c                the strict lower triangle is unaltered.
    ! c                if  info .ne. 0 , the factorization is not complete.
    ! c
    ! c        info    integer
    ! c                = 0  for normal return.
    ! c                = k  signals an error condition.  the leading minor
    ! c                     of order  k  is not positive definite.
    ! c
    ! c     linpack.  this version dated 08/14/78 .
    ! c     cleve moler, university of new mexico, argonne national lab.
    ! c      
    implicit none
    integer, intent(in) :: lda, n
    real (dp), intent(inout) :: a(lda,n)
    integer, intent(inout) :: info
    real (dp) :: t, s
    integer :: j,jm1,k

    ! c     begin block with ...exits to 40
    ! c
    ! c
    do j = 1, n
       info = j
       s = 0.0d0
       jm1 = j - 1
       if (jm1 < 1) go to 20
       do k = 1, jm1
          t = a(k,j) - ddot(k-1,a(1,k),1,a(1,j),1)
          t = t/a(k,k)
          a(k,j) = t
          s = s + t*t
       end do
20     continue
       s = a(j,j) - s
       ! c     ......exit          
       if (s <= 0.0d0) go to 40
       a(j,j) = sqrt(s)
    end do
    info = 0
40  continue
    return
  end subroutine dpofa

  function ddot(n,dx,incx,dy,incy)
    !
    !     forms the dot product of two vectors.
    !     uses unrolled loops for increments equal to one.
    !     jack dongarra, linpack, 3/11/78.
    !
    !     modified 12/3/93, array(1) declarations changed to array(*)
    !
    implicit none
    integer, intent(in) :: n, incx,incy
    real(dp), intent(in) :: dx(n*incx),dy(n*incy)
    integer :: i,ix,iy,m,mp1
    real(dp) :: ddot,dtemp

    ddot = 0.0d0
    dtemp = 0.0d0
    if(n.le.0)return
    if(incx.eq.1.and.incy.eq.1) go to 20

    ! c
    ! c        code for unequal increments or equal increments
    ! c          not equal to 1
    ! c
    ix = 1
    iy = 1
    if(incx.lt.0)ix = (-n+1)*incx + 1
    if(incy.lt.0)iy = (-n+1)*incy + 1
    do i = 1,n
       dtemp = dtemp + dx(ix)*dy(iy)
       ix = ix + incx
       iy = iy + incy
    end do
    ddot = dtemp
    return

    ! c
    ! c        code for both increments equal to 1
    ! c
    ! c
    ! c        clean-up loop
    ! c
20  m = mod(n,5)
    if( m .eq. 0 ) go to 40
    do i = 1,m
       dtemp = dtemp + dx(i)*dy(i)
    end do
    if( n .lt. 5 ) go to 60
40  mp1 = m + 1
    do i = mp1,n,5
       dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +&
            dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
    end do
60  ddot = dtemp
    return
  end function ddot

end module lpk
