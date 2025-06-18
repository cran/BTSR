module lib_utils
  !***************************************************************************
  !
  !  This module contains:
  !    - parameters
  !    - general purposes user defined variables
  !    - some special functions and subroutines
  !
  !***************************************************************************
  ! March 2025: added this module
  !***************************************************************************
  implicit none

  !---------------------------------------------------------------------
  ! The current models available
  !---------------------------------------------------------------------
  character(len=8), parameter :: model1 = "beta"
  character(len=8), parameter :: model2 = "gamma"
  character(len=8), parameter :: model3 = "kuma"
  character(len=8), parameter :: model4 = "matsu"
  character(len=8), parameter :: model5 = "ul"
  character(len=8), parameter :: model6 = "uw"
  ! array with the codes for each model
  character(len=8), parameter :: current_models(6) = &
       [model1, model2, model3, model4, model5, model6]
  ! array indicating if the model contains only nu
  logical, parameter :: onlymu(6) = &
       [.false., .false., .false., .true., .true., .false.]

  !---------------------------------------------------------------------
  ! some global parameters
  !---------------------------------------------------------------------
  integer,  parameter  :: dp = kind(1.d0)          ! double precision
  real(dp), parameter  :: pi = acos(-1.d0)         ! pi = 3.1415
  real(dp), parameter  :: infinity = 2*huge(1.0d0) ! Infinity
  real(dp), parameter  :: epsmch = epsilon(1.d0)   ! machine precision


  !-----------------------------------------------------------------------
  ! Interface - to use the same call for any allocation subroutine
  !-----------------------------------------------------------------------
  interface safe_allocate
     module procedure safe_allocateR1  ! vector 1:n
     module procedure safe_allocateR2  ! matrix 1:n1, 1:n2 or n1:m1, n2:m2
     module procedure safe_allocateI1  ! vector of integers
  end interface safe_allocate

contains

  !************************************************************************
  !
  !    Subroutines used to safe allocate vectors and matrices
  !
  !    Just to make sure that the program will not try to allocate
  !    a vector/matrix that is already allocated
  !
  !*************************************************************************
  subroutine safe_allocateR1(x, n1, n2)
    !*********************************************************************
    ! safe allocate real vectors.
    !
    ! Last revision: February, 2025
    ! - merged two subroutines and added the optional argument
    !*********************************************************************
    implicit none
    integer, intent(in) :: n1
    integer, intent(in), optional :: n2
    real(dp), allocatable, intent(inout) ::  x(:)

    !---------------------------------------------------------------------
    ! Deallocate if allocated
    !---------------------------------------------------------------------
    if (allocated(x)) deallocate(x)

    !---------------------------------------------------------------------
    ! Allocate with the specified bounds
    !---------------------------------------------------------------------
    if (present(n2)) then
       allocate(x(n1:n2))
    else
       allocate(x(n1))
    end if
    return
  end subroutine safe_allocateR1

  subroutine safe_allocateR2(array, n1, n2, n3, n4)
    !*********************************************************************
    ! safe allocate real matrices.
    !
    ! Last revision: February, 2025
    ! - merged two subroutines and added the optional argument
    !*********************************************************************
    implicit none
    real(dp), allocatable, intent(inout) :: array(:, :)
    integer, intent(in) :: n1, n2
    integer, intent(in), optional :: n3, n4

    !---------------------------------------------------------------------
    ! Deallocate if already allocated
    !---------------------------------------------------------------------
    if (allocated(array)) deallocate(array)

    !---------------------------------------------------------------------
    ! Allocate with the specified bounds
    !---------------------------------------------------------------------
    if (present(n3) .and. present(n4)) then
       allocate(array(n1:n2, n3:n4))
    else
       allocate(array(1:n1, 1:n2))
    end if
    return
  end subroutine safe_allocateR2

  subroutine safe_allocateI1(x, n1, n2)
    !---------------------------------------------------------------------
    ! safe allocate integer vectors
    !
    ! Last revision: February, 2025
    !  - added the optional paramter
    !---------------------------------------------------------------------
    implicit none
    integer, intent(in) :: n1
    integer, intent(in), optional :: n2
    integer, allocatable, intent(inout) ::  x(:)

    !---------------------------------------------------------------------
    ! Deallocate if allocated
    !---------------------------------------------------------------------
    if (allocated(x)) deallocate(x)

    !---------------------------------------------------------------------
    ! Allocate with the specified bounds
    !---------------------------------------------------------------------
    if (present(n2)) then
       allocate(x(n1:n2))
    else
       allocate(x(n1))
    end if
    return
  end subroutine safe_allocateI1

  !*************************************************************************
  !
  !                                 Special functions
  !
  !*************************************************************************
  function trigamma ( x, ifault )
    !*********************************************************************
    ! https://people.sc.fsu.edu/~jburkardt/f_src/asa121/asa121.f90
    !
    ! trigamma() calculates trigamma(x) = d^2 log(gamma(x)) / dx^2
    !
    !  Licensing:
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !    28 August 2021
    !
    !  Author:
    !    Original FORTRAN77 version by B.E. Schneider.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !    B.E. Schneider, Algorithm AS 121: Trigamma Function,
    !    Applied Statistics, Volume 27, Number 1, pages 97 - 99, 1978.
    !
    !  Input:
    !    real (dp) X, the argument of the trigamma function. 0 < X.
    !
    !  Output:
    !    integer IFAULT, error flag. 0, no error. 1, X <= 0.
    !    real (dp) TRIGAMMA, the value of the trigamma function.
    !
    !*********************************************************************
    implicit none
    real(dp), parameter :: a = 0.0001d+00
    real(dp), parameter :: b = 5.0d+00
    real(dp), parameter :: b2 =  0.1666666667d+00
    real(dp), parameter :: b4 =  - 0.03333333333d+00
    real(dp), parameter :: b6 =  0.02380952381d+00
    real(dp), parameter :: b8 =  - 0.03333333333d+00
    integer :: ifault
    real(dp) :: trigamma, x, y, z
    !---------------------------------------------------------------------
    !  Check the input.
    !---------------------------------------------------------------------
    if ( x <= 0.0d+00 ) then
       ifault = 1
       trigamma = 0.0d+00
       return
    end if

    ifault = 0
    z = x
    !---------------------------------------------------------------------
    !  Use small value approximation if X <= A.
    !---------------------------------------------------------------------
    if ( x <= a ) then
       trigamma = 1.0d+00 / x / x
       return
    end if
    !---------------------------------------------------------------------
    !  Increase argument to ( X + I ) >= B.
    !---------------------------------------------------------------------
    trigamma = 0.0d+00

    do while ( z < b )
       trigamma = trigamma + 1.0d+00 / z / z
       z = z + 1.0d+00
    end do
    !---------------------------------------------------------------------
    !  Apply asymptotic formula if argument is B or greater.
    !---------------------------------------------------------------------
    y = 1.0d+00 / z / z

    trigamma = trigamma + 0.5d+00 * y + (1.0d+00 + y * (b2 + y * (b4 + y * (b6 + y * b8)))) / z
    return
  end function trigamma

  function psi(xx) result(fn_val)
    !*********************************************************************
    !                 evaluation of the digamma function
    !                          -----------
    !     psi(xx) is assigned the value 0 when the digamma function cannot be computed.
    !     the main computation involves evaluation of rational chebyshev approximations
    !     published in math. comp. 27, 123 - 127(1973) by Cody, Strecok and Thacher.
    !*********************************************************************
    !     psi was written at Argonne National Laboratory for the funpack package of
    !     special function subroutines. psi was modified by A.H. Morris (nswc).
    !*********************************************************************
    implicit none
    real (dp), intent(in) :: xx
    real (dp)             :: fn_val
    real (dp) :: dx0 = 1.461632144968362341262659542325721325d0
    !---------------------------------------------------------------------
    !
    !     piov4 = pi / 4
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
    real (dp) :: p1(7) = [ .895385022981970d-02, .477762828042627d+01,  &
         .142441585084029d+03, .118645200713425d+04,  &
         .363351846806499d+04, .413810161269013d+04,  &
         .130560269827897d+04],   &
         q1(6) = [ .448452573429826d+02, .520752771467162d+03,  &
         .221000799247830d+04, .364127349079381d+04,  &
         .190831076596300d+04, .691091682714533d-05]
    !---------------------------------------------------------------------
    !
    !     coefficients for rational approximation of
    !     psi(x) - ln(x) + 1 / (2 * x),  x > 3.0
    !
    !---------------------------------------------------------------------
    real (dp) :: p2(4) = [-.212940445131011d+01,  -.701677227766759d+01,  &
         -.448616543918019d+01,  -.648157123766197d+00], &
         q2(4) = [.322703493791143d+02,  .892920700481861d+02,  &
         .546117738103215d+02,  .777788548522962d+01]
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
    !        xsmall = absolute argument below which pi * cotan(pi * x)
    !                 may be represented by 1 / x.
    !
    !---------------------------------------------------------------------
    xmax1 = Huge(1) ! the largest magnitude (huge)
    xmax1 = min(xmax1, 1.0d0 / epsmch) ! epsmch = the machine precision
    xsmall = 1.d-9
    !---------------------------------------------------------------------
    x = xx
    aug = 0.0d0
    if (x >= 0.5d0) go to 200
    !---------------------------------------------------------------------
    !     x .lt. 0.5,  use reflection formula
    !     psi(1 - x) = psi(x) + pi * cotan(pi * x)
    !---------------------------------------------------------------------
    if (abs(x) > xsmall) go to 100
    if (x == 0.0d0) go to 400
    !---------------------------------------------------------------------
    !     0 .lt. abs(x) .le. xsmall.  use 1 / x as a substitute
    !     for  pi * cotan(pi * x)
    !---------------------------------------------------------------------
    aug =  - 1.0d0 / x
    go to 150
    !---------------------------------------------------------------------
    !     reduction of argument for cotan
    !---------------------------------------------------------------------
100 w = - x
    sgn = piov4
    if (w > 0.0d0) go to 120
    w = - w
    sgn =  - sgn
    !---------------------------------------------------------------------
    !     make an error exit if x .le.  - xmax1
    !---------------------------------------------------------------------
120 if (w >= xmax1) go to 400
    nq = int(w)
    w = w - nq
    nq = int(w * 4.0d0)
    w = 4.0d0 * (w - nq * .25d0)
    !---------------------------------------------------------------------
    !     w is now related to the fractional part of  4.0 * x.
    !     adjust argument to correspond to values in first
    !     quadrant and determine sign
    !---------------------------------------------------------------------
    n = nq / 2
    if ((n + n)  /= nq) w = 1.0d0 - w
    z = piov4 * w
    m = n / 2
    if ((m + m)  /= n) sgn = - sgn
    !---------------------------------------------------------------------
    !     determine final value for   - pi * cotan(pi * x)
    !---------------------------------------------------------------------
    n = (nq + 1) / 2
    m = n / 2
    m = m + m
    if (m  /= n) go to 140
    !---------------------------------------------------------------------
    !     check for singularity
    !---------------------------------------------------------------------
    if (z == 0.0d0) go to 400
    !---------------------------------------------------------------------
    !     use cos / sin as a substitute for cotan, and
    !     sin / cos as a substitute for tan
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
       upper = (upper + p1(i + 1)) * x
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
       upper = (upper + p2(i + 1)) * w
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

  function digamma(x) result(fn)
    !*********************************************************************
    ! digamma function
    ! Calculates the value of psi(x)
    !
    ! Implemented by Taiane S. Prass
    ! March, 2018
    !*********************************************************************
    implicit none
    real(dp), intent(in) :: x
    real (dp) :: fn
    fn = psi(x)
  end function digamma

  elemental logical function is_finite(x)
    !*********************************************************************
    ! This function checks for nan and non finite values
    !
    ! Implemented by Taiane S. Prass
    ! February, 2025
    !*********************************************************************
    implicit none
    real(dp), intent(in) :: x
    is_finite = (x < huge(x)) .and. (x > -huge(x)) .and. (x == x)
    return
  end function is_finite

end module lib_utils
