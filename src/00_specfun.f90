module specfun
  !----------------------------------------------
  !
  !  Special functions needed by subroutines
  !
  !  Codes were colected from diferent modules.
  !  See each function for details.
  !
  !----------------------------------------------

  ! double precision:    
  private :: dp
  integer, parameter :: dp = kind(1.d0)   

contains   

  function trigamma ( x, ifault )

    !*****************************************************************************80
    !
    ! https://people.sc.fsu.edu/~jburkardt/f_src/asa121/asa121.f90
    !
    ! trigamma() calculates trigamma(x) = d^2 log(gamma(x)) / dx^2
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    28 August 2021
    !
    !  Author:
    !
    !    Original FORTRAN77 version by B.E. Schneider.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    B.E. Schneider,
    !    Algorithm AS 121:
    !    Trigamma Function,
    !    Applied Statistics, 
    !    Volume 27, Number 1, pages 97-99, 1978.
    !
    !  Input:
    !
    !    real (dp) X, the argument of the trigamma function.
    !    0 < X.
    !
    !  Output:
    !
    !    integer IFAULT, error flag.
    !    0, no error.
    !    1, X <= 0.
    !
    !    real (dp) TRIGAMMA, the value of the trigamma function.
    ! 
    ! integer, parameter :: rk = kind ( 1.0D+00 )
    !
    implicit none
    real(dp), parameter :: a = 0.0001D+00
    real(dp), parameter :: b = 5.0D+00
    real(dp), parameter :: b2 =  0.1666666667D+00
    real(dp), parameter :: b4 = -0.03333333333D+00
    real(dp), parameter :: b6 =  0.02380952381D+00
    real(dp), parameter :: b8 = -0.03333333333D+00
    integer :: ifault
    real(dp) :: trigamma
    real(dp) :: x
    real(dp) :: y
    real(dp) :: z
    !
    !  Check the input.
    !
    if ( x <= 0.0D+00 ) then
       ifault = 1
       trigamma = 0.0D+00
       return
    end if

    ifault = 0
    z = x
    !
    !  Use small value approximation if X <= A.
    !
    if ( x <= a ) then
       trigamma = 1.0D+00 / x / x
       return
    end if
    !
    !  Increase argument to ( X + I ) >= B.
    !
    trigamma = 0.0D+00

    do while ( z < b )
       trigamma = trigamma + 1.0D+00 / z / z
       z = z + 1.0D+00
    end do
    !
    !  Apply asymptotic formula if argument is B or greater.
    !
    y = 1.0D+00 / z / z

    trigamma = trigamma + 0.5D+00 * &
         y + ( 1.0D+00 &
         + y * ( b2  &
         + y * ( b4  &
         + y * ( b6  &
         + y *   b8 )))) / z

    return
  end function trigamma

  function lngamma(z) result(lanczos)
    !
    ! https://jblevins.org/mirror/amiller/lanczos.f90
    !
    !  Uses Lanczos-type approximation to ln(gamma) for z > 0.
    !
    !  Reference:
    !       Lanczos, C. 'A precision approximation of the gamma
    !               function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
    !  Accuracy: About 14 significant digits except for small regions
    !            in the vicinity of 1 and 2.
    !  Programmer: Alan Miller
    !              1 Creswick Street, Brighton, Vic. 3187, Australia
    !  e-mail: amiller @ bigpond.net.au
    !  Latest revision - 14 October 1996
    implicit none
    real(dp), intent(in) :: z
    real(dp) :: lanczos
    ! local variables
    real(dp)  :: a(9) = (/ 0.9999999999995183d0, 676.5203681218835d0, &
         -1259.139216722289d0, 771.3234287757674d0, &
         -176.6150291498386d0, 12.50734324009056d0, &
         -0.1385710331296526d0, 0.9934937113930748d-05, &
         0.1659470187408462d-06 /), zero = 0.d0,   &
         one = 1.d0, lnsqrt2pi = 0.9189385332046727d0, &
         half = 0.5d0, sixpt5 = 6.5d0, seven = 7.d0, tmp
    integer   :: j

    lanczos = zero
    if (z <= zero) then
       !write(*, *) 'error: zero or -ve argument for lngamma'
       return
    end if
    tmp = z + seven
    do j = 9, 2, -1
       lanczos = lanczos + a(j)/tmp
       tmp = tmp - one
    end do
    lanczos = lanczos + a(1)
    lanczos = log(lanczos) + lnsqrt2pi - (z + sixpt5) + (z - half)*log(z + sixpt5)
    return
  end function lngamma

  function alnrel(a) result(fn_val)
    !-----------------------------------------------------------------------
    !            evaluation of the function ln(1 + a)
    !-----------------------------------------------------------------------
    implicit none
    real (dp), intent(in) :: a
    real (dp)             :: fn_val

    ! local variables
    real (dp) :: p1 = -.129418923021993d+01, p2 = .405303492862024d+00,  &
         p3 = -.178874546012214d-01, q1 = -.162752256355323d+01, &
         q2 = .747811014037616d+00, q3 = -.845104217945565d-01,  &
         t, t2, w, x, zero = 0.d0, half = 0.5d0, one = 1.d0, two = 2.d0
    !--------------------------
    if (abs(a) <= 0.375d0) then
       t = a/(a + two)
       t2 = t*t
       w = (((p3*t2 + p2)*t2 + p1)*t2 + one)/ (((q3*t2 + q2)*t2 + q1)*t2 + one)
       fn_val = two*t*w
    else
       x = one + a
       if (a < zero) x = (a + half) + half
       fn_val = log(x)
    end if

    return
  end function alnrel

  function algdiv (a, b) result(fn_val)
    !-----------------------------------------------------------------------
    !
    !     computation of ln(gamma(b)/gamma(a+b)) when b >= 8
    !
    !                         --------
    !
    !     in this algorithm, del(x) is the function defined by
    !     ln(gamma(x)) = (x - 0.5)*ln(x) - x + 0.5*ln(2*pi) + del(x).
    !
    !-----------------------------------------------------------------------
    implicit none
    real (dp), intent(in) :: a, b
    real (dp)             :: fn_val

    ! external   alnrel
    real (dp) :: c, d, h, s11, s3, s5, s7, s9, t, u, v, w, x, x2
    real (dp) :: c0 = .833333333333333d-01, c1 = -.277777777760991d-02, &
         c2 = .793650666825390d-03, c3 = -.595202931351870d-03, &
         c4 = .837308034031215d-03, c5 = -.165322962780713d-02
    !------------------------
    if (a > b) then
       h = b/a
       c = 1.0d0/(1.0d0 + h)
       x = h/(1.0d0 + h)
       d = a + (b - 0.5d0)
    else
       h = a/b
       c = h/(1.0d0 + h)
       x = 1.0d0/(1.0d0 + h)
       d = b + (a - 0.5d0)
    end if
    !
    !                set sn = (1 - x**n)/(1 - x)
    !
    x2 = x*x
    s3 = 1.0d0 + (x + x2)
    s5 = 1.0d0 + (x + x2*s3)
    s7 = 1.0d0 + (x + x2*s5)
    s9 = 1.0d0 + (x + x2*s7)
    s11 = 1.0d0 + (x + x2*s9)
    !
    !                set w = del(b) - del(a + b)
    !
    t = (1.0d0/b)**2
    w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3)*t + c0
    w = w*(c/b)
    !
    !                    combine the results
    !
    u = d*alnrel(a/b)
    v = a*(log(b) - 1.0d0)
    if (u > v) then
       fn_val = (w - v) - u
    else
       fn_val = (w - u) - v
    end if

    return
  end function algdiv

  function gsumln (a, b) result(fn_val)
    !-----------------------------------------------------------------------
    !          evaluation of the function ln(gamma(a + b))
    !          for 1 .le. a .le. 2  and  1 .le. b .le. 2
    !-----------------------------------------------------------------------
    implicit none
    real (dp), intent(in) :: a, b
    real (dp)             :: fn_val

    real (dp) :: x
    x = a + b - 2.d0
    if (x > 0.25d0) go to 10
    fn_val = gamln1(1.0d0 + x)
    return
10  if (x > 1.25d0) go to 20
    fn_val = gamln1(x) + alnrel(x)
    return
20  fn_val = gamln1(x - 1.0d0) + log(x*(1.0d0 + x))
    return
  end function gsumln

  function bcorr (a0, b0) result(fn_val)
    !-----------------------------------------------------------------------
    !
    !     evaluation of  del(a0) + del(b0) - del(a0 + b0)  where
    !     ln(gamma(a)) = (a - 0.5)*ln(a) - a + 0.5*ln(2*pi) + del(a).
    !     it is assumed that a0 .ge. 8 and b0 .ge. 8.
    !
    !-----------------------------------------------------------------------
    implicit none
    real (dp), intent(in) :: a0, b0
    real (dp)             :: fn_val

    real (dp) :: a, b, c, h, s11, s3, s5, s7, s9, t, w, x, x2,         &
         c0 = .833333333333333d-01, c1 = -.277777777760991d-02,  &
         c2 = .793650666825390d-03, c3 = -.595202931351870d-03,  &
         c4 = .837308034031215d-03, c5 = -.165322962780713d-02
    !------------------------
    a = min(a0, b0)
    b = max(a0, b0)

    h = a/b
    c = h/(1.0d0 + h)
    x = 1.0d0/(1.0d0 + h)
    x2 = x*x
    !
    !                set sn = (1 - x**n)/(1 - x)
    !
    s3 = 1.0d0 + (x + x2)
    s5 = 1.0d0 + (x + x2*s3)
    s7 = 1.0d0 + (x + x2*s5)
    s9 = 1.0d0 + (x + x2*s7)
    s11 = 1.0d0 + (x + x2*s9)
    !
    !                set w = del(b) - del(a + b)
    !
    t = (1.0d0/b)**2
    w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3)*t + c0
    w = w*(c/b)
    !
    !                   compute  del(a) + w
    !
    t = (1.0d0/a)**2
    fn_val = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0)/a + w
    return
  end function bcorr

  function betaln (a0, b0) result(fn_val)
    !-----------------------------------------------------------------------
    !     evaluation of the logarithm of the beta function
    !-----------------------------------------------------------------------
    !     e = 0.5*ln(2*pi)
    !--------------------------
    implicit none
    real (dp), intent(in) :: a0, b0
    real (dp)             :: fn_val

    real (dp) :: a, b, c, e = .918938533204673d0, h, u, v, w, z
    integer   :: i, n
    !--------------------------
    a = min(a0,b0)
    b = max(a0,b0)
    if (a >= 8.0d0) go to 60
    if (a >= 1.0d0) go to 20
    !-----------------------------------------------------------------------
    !                   procedure when a .lt. 1
    !-----------------------------------------------------------------------
    if (b >= 8.0d0) go to 10
    fn_val = gamln(a) + (gamln(b) - gamln(a + b))
    return
10  fn_val = gamln(a) + algdiv(a,b)
    return
    !-----------------------------------------------------------------------
    !                procedure when 1 .le. a .lt. 8
    !-----------------------------------------------------------------------
20  if (a > 2.0d0) go to 30
    if (b > 2.0d0) go to 21
    fn_val = gamln(a) + gamln(b) - gsumln(a,b)
    return
21  w = 0.0d0
    if (b < 8.0d0) go to 40
    fn_val = gamln(a) + algdiv(a,b)
    return
    !
    !                reduction of a when b .le. 1000
    !
30  if (b > 1000.0d0) go to 50
    n = int(a - 1.0d0)
    w = 1.0d0
    do i = 1, n
       a = a - 1.0d0
       h = a/b
       w = w * (h/(1.0d0 + h))
    end do
    w = log(w)
    if (b < 8.0d0) go to 40
    fn_val = w + gamln(a) + algdiv(a,b)
    return
    !
    !                 reduction of b when b .lt. 8
    !
40  n = int(b - 1.0d0)
    z = 1.0d0
    do i = 1,n
       b = b - 1.0d0
       z = z * (b/(a + b))
    end do
    fn_val = w + log(z) + (gamln(a) + (gamln(b) - gsumln(a,b)))
    return
    !
    !                reduction of a when b .gt. 1000
    !
50  n = int(a - 1.0d0)
    w = 1.0d0
    do i = 1,n
       a = a - 1.0d0
       w = w * (a/(1.0d0 + a/b))
    end do
    fn_val = (log(w) - n*log(b)) + (gamln(a) + algdiv(a,b))
    return
    !-----------------------------------------------------------------------
    !                   procedure when a .ge. 8
    !-----------------------------------------------------------------------
60  w = bcorr(a,b)
    h = a/b
    c = h/(1.0d0 + h)
    u = -(a - 0.5d0)*log(c)
    v = b*alnrel(h)
    if (u <= v) go to 61
    fn_val = (((-0.5d0*log(b) + e) + w) - v) - u
    return
61  fn_val = (((-0.5d0*log(b) + e) + w) - u) - v
    return
  end function betaln

  function rlog1(x) result(fn_val)
    !-----------------------------------------------------------------------
    !             evaluation of the function x - ln(1 + x)
    !-----------------------------------------------------------------------
    !     a = rlog (0.7)
    !     b = rlog (4/3)
    !------------------------
    implicit none
    real (dp), intent(in) :: x
    real (dp)             :: fn_val

    real (dp) ::   a = .566749439387324d-01, b = .456512608815524d-01,  &
         p0 = .333333333333333d+00, p1 = -.224696413112536d+00, &
         p2 = .620886815375787d-02, q1 = -.127408923933623d+01, &
         q2 = .354508718369557d+00, r, t, u, up2, w, w1
    !------------------------
    if (x < -0.39d0 .or. x > 0.57d0) go to 100
    if (x < -0.18d0) go to 10
    if (x >  0.18d0) go to 20
    !
    !                 argument reduction
    !
    u = x
    up2 = u + 2.0d0
    w1 = 0.0d0
    go to 30

10  u = (x + 0.3d0)/0.7d0
    up2 = u + 2.0d0
    w1 = a - u*0.3d0
    go to 30

20  t = 0.75d0*x
    u = t - 0.25d0
    up2 = t + 1.75d0
    w1 = b + u/3.0d0
    !
    !                  series expansion
    !
30  r = u/up2
    t = r*r
    w = ((p2*t + p1)*t + p0)/((q2*t + q1)*t + 1.0d0)
    fn_val = r*(u - 2.0d0*t*w) + w1
    return

100 w = (x + 0.5d0) + 0.5d0
    fn_val = x - log(w)
    return
  end function rlog1

  function gamln1 (a) result(fn_val)
    !-----------------------------------------------------------------------
    !     evaluation of ln(gamma(1 + a)) for -0.2 .le. a .le. 1.25
    !-----------------------------------------------------------------------
    implicit none
    real (dp), intent(in) :: a
    real (dp)             :: fn_val

    real (dp) :: w, x, &
         p0 =  .577215664901533d+00, p1 =  .844203922187225d+00,  &
         p2 = -.168860593646662d+00, p3 = -.780427615533591d+00,  &
         p4 = -.402055799310489d+00, p5 = -.673562214325671d-01,  &
         p6 = -.271935708322958d-02,   &
         q1 =  .288743195473681d+01, q2 =  .312755088914843d+01,  &
         q3 =  .156875193295039d+01, q4 =  .361951990101499d+00,  &
         q5 =  .325038868253937d-01, q6 =  .667465618796164d-03,  &
         r0 = .422784335098467d+00,  r1 = .848044614534529d+00,  &
         r2 = .565221050691933d+00,  r3 = .156513060486551d+00,  &
         r4 = .170502484022650d-01,  r5 = .497958207639485d-03,  &
         s1 = .124313399877507d+01,  s2 = .548042109832463d+00,  &
         s3 = .101552187439830d+00,  s4 = .713309612391000d-02,  &
         s5 = .116165475989616d-03
    !----------------------
    if (a >= 0.6d0) go to 10
    w = ((((((p6*a + p5)*a + p4)*a + p3)*a + p2)*a + p1)*a + p0)/  &
         ((((((q6*a + q5)*a + q4)*a + q3)*a + q2)*a + q1)*a + 1.0d0)
    fn_val = -a*w
    return

10  x = (a - 0.5d0) - 0.5d0
    w = (((((r5*x + r4)*x + r3)*x + r2)*x + r1)*x + r0)/  &
         (((((s5*x + s4)*x + s3)*x + s2)*x + s1)*x + 1.0d0)
    fn_val = x*w
    return
  end function gamln1

  function gamln (a) result(fn_val)
    !-----------------------------------------------------------------------
    !            evaluation of ln(gamma(a)) for positive a
    !-----------------------------------------------------------------------
    !     written by alfred h. morris
    !          naval surface warfare center
    !          dahlgren, virginia
    !--------------------------
    !     d = 0.5*(ln(2*pi) - 1)
    !--------------------------
    implicit none
    real (dp), intent(in) :: a
    real (dp)             :: fn_val

    real (dp) :: c0 = .833333333333333d-01, c1 = -.277777777760991d-02,  &
         c2 = .793650666825390d-03, c3 = -.595202931351870d-03,  &
         c4 = .837308034031215d-03, c5 = -.165322962780713d-02,  &
         d = .418938533204673d0, t, w
    integer   :: i, n
    !--------------------------
    if (a > 0.8d0) go to 10
    fn_val = gamln1(a) - log(a)
    return
10  if (a > 2.25d0) go to 20
    t = (a - 0.5d0) - 0.5d0
    fn_val = gamln1(t)
    return

20  if (a >= 10.0d0) go to 30
    n = int(a - 1.25d0)
    t = a
    w = 1.0d0
    do i = 1, n
       t = t - 1.0d0
       w = t*w
    end do
    fn_val = gamln1(t - 1.0d0) + log(w)
    return

30  t = (1.0d0/a)**2
    w = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0)/a
    fn_val = (d + w) + (a - 0.5d0)*(log(a) - 1.0d0)

    return
  end function gamln

  function gam1(a) result(fn_val)
    !-----------------------------------------------------------------------
    !     computation of 1/gamma(a+1) - 1  for -0.5 <= a <= 1.5
    !-----------------------------------------------------------------------
    implicit none
    real (dp), intent(in) :: a
    real (dp)             :: fn_val

    real (dp) :: p(7) = (/  .577215664901533d+00, -.409078193005776d+00, &
         -.230975380857675d+00,  .597275330452234d-01, &
         .766968181649490d-02, -.514889771323592d-02, &
         .589597428611429d-03 /),   &
         q(5) = (/  .100000000000000d+01,  .427569613095214d+00, &
         .158451672430138d+00,  .261132021441447d-01, &
         .423244297896961d-02 /),   &
         r(9) = (/ -.422784335098468d+00, -.771330383816272d+00, &
         -.244757765222226d+00,  .118378989872749d+00, &
         .930357293360349d-03, -.118290993445146d-01, &
         .223047661158249d-02,  .266505979058923d-03, &
         -.132674909766242d-03 /)
    !------------------------
    real (dp) :: bot, d, s1 = .273076135303957d+00, s2 = .559398236957378d-01,  &
         t, top, w
    !------------------------
    t = a
    d = a - 0.5d0
    if (d > 0.0d0) t = d - 0.5d0

    if (t > 0.d0) then
       top = (((((p(7)*t + p(6))*t + p(5))*t + p(4))*t + p(3))*t + p(2))*t + p(1)
       bot = (((q(5)*t + q(4))*t + q(3))*t + q(2))*t + 1.0d0
       w = top/bot
       if (d > 0.0d0) then
          fn_val = (t/a)*((w - 0.5d0) - 0.5d0)
       else
          fn_val = a*w
       end if
    else if (t < 0.d0) then
       top = (((((((r(9)*t + r(8))*t + r(7))*t + r(6))*t + r(5))*t  &
            + r(4))*t + r(3))*t + r(2))*t + r(1)
       bot = (s2*t + s1)*t + 1.0d0
       w = top/bot
       if (d > 0.0d0) then
          fn_val = t*w/a
       else
          fn_val = a*((w + 0.5d0) + 0.5d0)
       end if
    else
       fn_val = 0.0d0
    end if

    return
  end function gam1

  function brcomp (a, b, x, y) result(fn_val)
    !-----------------------------------------------------------------------
    !               evaluation of x**a*y**b/beta(a,b)
    !-----------------------------------------------------------------------
    implicit none
    real (dp), intent(in) :: a, b, x, y
    real (dp)             :: fn_val

    real (dp) :: lambda, lnx, lny
    !-----------------
    !     const = 1/sqrt(2*pi)
    !-----------------
    real (dp) :: apb, a0, b0, c, const = .398942280401433d0, e, h, t, u,  &
         v, x0, y0, z
    integer   :: i, n

    fn_val = 0.0d0
    if (x == 0.0d0 .or. y == 0.0d0) return
    a0 = min(a,b)
    if (a0 >= 8.0d0) go to 100

    if (x > 0.375d0) go to 10
    lnx = log(x)
    lny = alnrel(-x)
    go to 20
10  if (y > 0.375d0) go to 11
    lnx = alnrel(-y)
    lny = log(y)
    go to 20
11  lnx = log(x)
    lny = log(y)

20  z = a*lnx + b*lny
    if (a0 < 1.0d0) go to 30
    z = z - betaln(a,b)
    fn_val = exp(z)
    return
    !-----------------------------------------------------------------------
    !              procedure for a .lt. 1 or b .lt. 1
    !-----------------------------------------------------------------------
30  b0 = max(a,b)
    if (b0 >= 8.0d0) go to 80
    if (b0 > 1.0d0) go to 60
    !
    !                   algorithm for b0 .le. 1
    !
    fn_val = exp(z)
    if (fn_val == 0.0d0) return

    apb = a + b
    if (apb > 1.0d0) go to 40
    z = 1.0d0 + gam1(apb)
    go to 50
40  u = a + b - 1.d0
    z = (1.0d0 + gam1(u))/apb

50  c = (1.0d0 + gam1(a))*(1.0d0 + gam1(b))/z
    fn_val = fn_val*(a0*c)/(1.0d0 + a0/b0)
    return
    !
    !                algorithm for 1 .lt. b0 .lt. 8
    !
60  u = gamln1(a0)
    n = int(b0 - 1.0d0)
    if (n < 1) go to 70
    c = 1.0d0
    do i = 1, n
       b0 = b0 - 1.0d0
       c = c*(b0/(a0 + b0))
    end do
    u = log(c) + u

70  z = z - u
    b0 = b0 - 1.0d0
    apb = a0 + b0
    if (apb > 1.0d0) go to 71
    t = 1.0d0 + gam1(apb)
    go to 72
71  u = a0 + b0 - 1.d0
    t = (1.0d0 + gam1(u))/apb
72  fn_val = a0*exp(z)*(1.0d0 + gam1(b0))/t
    return
    !
    !                   algorithm for b0 .ge. 8
    !
80  u = gamln1(a0) + algdiv(a0,b0)
    fn_val = a0*exp(z - u)
    return
    !-----------------------------------------------------------------------
    !              procedure for a .ge. 8 and b .ge. 8
    !-----------------------------------------------------------------------
100 if (a > b) go to 101
    h = a/b
    x0 = h/(1.0d0 + h)
    y0 = 1.0d0/(1.0d0 + h)
    lambda = a - (a + b)*x
    go to 110
101 h = b/a
    x0 = 1.0d0/(1.0d0 + h)
    y0 = h/(1.0d0 + h)
    lambda = (a + b)*y - b

110 e = -lambda/a
    if (abs(e) > 0.6d0) go to 111
    u = rlog1(e)
    go to 120
111 u = e - log(x/x0)

120 e = lambda/b
    if (abs(e) > 0.6d0) go to 121
    v = rlog1(e)
    go to 130
121 v = e - log(y/y0)

130 z = exp(-(a*u + b*v))
    fn_val = const*sqrt(b*x0)*z*exp(-bcorr(a,b))
    return
  end function brcomp

  function ipmpar (i) result(fn_val)
    !-----------------------------------------------------------------------
    !
    !     ipmpar provides the integer machine constants for the computer
    !     that is used. it is assumed that the argument i is an integer
    !     having one of the values 1-10. ipmpar(i) has the value ...
    !
    !  integers.
    !
    !     assume integers are represented in the n-digit, base-a form
    !
    !               sign ( x(n-1)*a**(n-1) + ... + x(1)*a + x(0) )
    !
    !               where 0 .le. x(i) .lt. a for i=0,...,n-1.
    !
    !     ipmpar(1) = a, the base (radix).
    !
    !     ipmpar(2) = n, the number of base-a digits (digits).
    !
    !     ipmpar(3) = a**n - 1, the largest magnitude (huge).
    !
    !  floating-point numbers.
    !
    !     it is assumed that the single and real floating
    !     point arithmetics have the same base, say b, and that the
    !     nonzero numbers are represented in the form
    !
    !               sign (b**e) * (x(1)/b + ... + x(m)/b**m)
    !
    !               where x(i) = 0,1,...,b-1 for i=1,...,m,
    !               x(1) .ge. 1, and emin .le. e .le. emax.
    !
    !     ipmpar(4) = b, the base.
    !
    !  single-precision
    !
    !     ipmpar(5) = m, the number of base-b digits.
    !
    !     ipmpar(6) = emin, the smallest exponent e.
    !
    !     ipmpar(7) = emax, the largest exponent e.
    !
    !  double-precision
    !
    !     ipmpar(8) = m, the number of base-b digits.
    !
    !     ipmpar(9) = emin, the smallest exponent e.
    !
    !     ipmpar(10) = emax, the largest exponent e.
    !
    !-----------------------------------------------------------------------

    integer, intent(in) :: i
    integer             :: fn_val

    select case(i)
    case(1)
       fn_val = radix(i)
    case(2)
       fn_val = digits(i)
    case(3)
       fn_val = huge(i)
    case(4)
       fn_val = radix(1.0)
    case(5)
       fn_val = digits(1.0)
    case(6)
       fn_val = minexponent(1.0)
    case(7)
       fn_val = maxexponent(1.0)
    case(8)
       fn_val = digits(1.0d0)
    case(9)
       fn_val = minexponent(1.0d0)
    case(10)
       fn_val = maxexponent(1.0d0)
    case default
       fn_val = 0.d0
    end select

    return
  end function ipmpar

  function dpmpar (i) result(fn_val)
    !-----------------------------------------------------------------------
    !
    !     dpmpar provides the real machine constants for
    !     the computer being used. it is assumed that the argument
    !     i is an integer having one of the values 1, 2, or 3. if the
    !     real arithmetic being used has m base b digits and
    !     its smallest and largest exponents are emin and emax, then
    !
    !        dpmpar(1) = b**(1 - m), the machine precision,
    !
    !        dpmpar(2) = b**(emin - 1), the smallest magnitude,
    !
    !        dpmpar(3) = b**emax*(1 - b**(-m)), the largest magnitude.
    !-----------------------------------------------------------------------

    integer, intent(in) :: i
    real (dp)           :: fn_val

    ! local variable
    real (dp)    :: one = 1.d0

    select case (i)
    case (1)
       fn_val = epsilon(one)
    case (2)
       fn_val = tiny(one)
    case (3)
       fn_val = huge(one)
    case default
       fn_val = 0.d0
    end select

    return
  end function dpmpar

  function psi(xx) result(fn_val)
    !---------------------------------------------------------------------
    !
    !                 evaluation of the digamma function
    !
    !                           -----------
    !
    !     psi(xx) is assigned the value 0 when the digamma function cannot
    !     be computed.
    !
    !     the main computation involves evaluation of rational chebyshev
    !     approximations published in math. comp. 27, 123-127(1973) by
    !     cody, strecok and thacher.
    !
    !---------------------------------------------------------------------
    !     psi was written at argonne national laboratory for the funpack
    !     package of special function subroutines. psi was modified by
    !     a.h. morris (nswc).
    !---------------------------------------------------------------------
    implicit none
    real (dp), intent(in) :: xx
    real (dp)             :: fn_val

    real (dp) :: dx0 = 1.461632144968362341262659542325721325d0
    !---------------------------------------------------------------------
    !
    !     piov4 = pi/4
    !     dx0 = zero of psi to extended precision
    !
    !---------------------------------------------------------------------
    real (dp) :: aug, den, piov4 = .785398163397448d0, sgn, upper,  &
         w, x, xmax1, xmx0, xsmall, z
    integer   :: i, m, n, nq
    !---------------------------------------------------------------------
    !
    !     coefficients for rational approximation of
    !     psi(x) / (x - x0),  0.5 <= x <= 3.0
    !
    !---------------------------------------------------------------------
    real (dp) :: p1(7) = (/ .895385022981970d-02, .477762828042627d+01,  &
         .142441585084029d+03, .118645200713425d+04,  &
         .363351846806499d+04, .413810161269013d+04,  &
         .130560269827897d+04 /),   &
         q1(6) = (/ .448452573429826d+02, .520752771467162d+03,  &
         .221000799247830d+04, .364127349079381d+04,  &
         .190831076596300d+04, .691091682714533d-05 /)
    !---------------------------------------------------------------------
    !
    !     coefficients for rational approximation of
    !     psi(x) - ln(x) + 1 / (2*x),  x > 3.0
    !
    !---------------------------------------------------------------------
    real (dp) :: p2(4) = (/ -.212940445131011d+01, -.701677227766759d+01,  &
         -.448616543918019d+01, -.648157123766197d+00 /), &
         q2(4) = (/  .322703493791143d+02,  .892920700481861d+02,  &
         .546117738103215d+02,  .777788548522962d+01 /)
    !---------------------------------------------------------------------
    !
    !     machine dependent constants ...
    !
    !        xmax1  = the smallest positive floating point constant
    !                 with entirely integer representation.  also used
    !                 as negative of lower bound on acceptable negative
    !                 arguments and as the positive argument beyond which
    !                 psi may be represented as alog(x).
    !
    !        xsmall = absolute argument below which pi*cotan(pi*x)
    !                 may be represented by 1/x.
    !
    !---------------------------------------------------------------------
    xmax1 = ipmpar(3)
    xmax1 = min(xmax1, 1.0d0/dpmpar(1))
    xsmall = 1.d-9
    !---------------------------------------------------------------------
    x = xx
    aug = 0.0d0
    if (x >= 0.5d0) go to 200
    !---------------------------------------------------------------------
    !     x .lt. 0.5,  use reflection formula
    !     psi(1-x) = psi(x) + pi * cotan(pi*x)
    !---------------------------------------------------------------------
    if (abs(x) > xsmall) go to 100
    if (x == 0.0d0) go to 400
    !---------------------------------------------------------------------
    !     0 .lt. abs(x) .le. xsmall.  use 1/x as a substitute
    !     for  pi*cotan(pi*x)
    !---------------------------------------------------------------------
    aug = -1.0d0 / x
    go to 150
    !---------------------------------------------------------------------
    !     reduction of argument for cotan
    !---------------------------------------------------------------------
100 w = - x
    sgn = piov4
    if (w > 0.0d0) go to 120
    w = - w
    sgn = -sgn
    !---------------------------------------------------------------------
    !     make an error exit if x .le. -xmax1
    !---------------------------------------------------------------------
120 if (w >= xmax1) go to 400
    nq = int(w)
    w = w - nq
    nq = int(w*4.0d0)
    w = 4.0d0 * (w - nq * .25d0)
    !---------------------------------------------------------------------
    !     w is now related to the fractional part of  4.0 * x.
    !     adjust argument to correspond to values in first
    !     quadrant and determine sign
    !---------------------------------------------------------------------
    n = nq / 2
    if ((n+n) /= nq) w = 1.0d0 - w
    z = piov4 * w
    m = n / 2
    if ((m+m) /= n) sgn = - sgn
    !---------------------------------------------------------------------
    !     determine final value for  -pi*cotan(pi*x)
    !---------------------------------------------------------------------
    n = (nq + 1) / 2
    m = n / 2
    m = m + m
    if (m /= n) go to 140
    !---------------------------------------------------------------------
    !     check for singularity
    !---------------------------------------------------------------------
    if (z == 0.0d0) go to 400
    !---------------------------------------------------------------------
    !     use cos/sin as a substitute for cotan, and
    !     sin/cos as a substitute for tan
    !---------------------------------------------------------------------
    aug = sgn * ((cos(z) / sin(z)) * 4.0d0)
    go to 150
140 aug = sgn * ((sin(z) / cos(z)) * 4.0d0)
150 x = 1.0d0 - x
200 if (x > 3.0d0) go to 300
    !---------------------------------------------------------------------
    !     0.5 .le. x .le. 3.0
    !---------------------------------------------------------------------
    den = x
    upper = p1(1) * x

    do i = 1, 5
       den = (den + q1(i)) * x
       upper = (upper + p1(i+1)) * x
    end do

    den = (upper + p1(7)) / (den + q1(6))
    xmx0 = x - dx0
    fn_val = den * xmx0 + aug
    return
    !---------------------------------------------------------------------
    !     if x .ge. xmax1, psi = ln(x)
    !---------------------------------------------------------------------
300 if (x >= xmax1) go to 350
    !---------------------------------------------------------------------
    !     3.0 .lt. x .lt. xmax1
    !---------------------------------------------------------------------
    w = 1.0d0 / (x * x)
    den = w
    upper = p2(1) * w

    do i = 1, 3
       den = (den + q2(i)) * w
       upper = (upper + p2(i+1)) * w
    end do

    aug = upper / (den + q2(4)) - 0.5d0 / x + aug
350 fn_val = aug + log(x)
    return
    !---------------------------------------------------------------------
    !     error return
    !---------------------------------------------------------------------
400 fn_val = 0.0d0
    return
  end function psi

end module specfun

