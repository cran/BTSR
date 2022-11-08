module RNG_mod
  !----------------------------------------------------------
  !
  !   This module contains subroutines related to
  ! random number generation.
  !
  !   The subroutines were colected from different sources.
  ! References are provided with each subroutine. 
  !
  !----------------------------------------------------------

  use main_mod, only : argsDist
  use specfun
  implicit none
  private
  public :: rng_t, rng_seed, rng_uniform
  public :: rbeta, dbeta
  public :: rkuma, dkuma
  public :: rgamma, d_gamma
  public :: ruw, duw

  integer, parameter :: dp = kind(1.d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383280d0
  real(dp), parameter :: M_2PI = 6.283185307179586476925286766559d0   !/* 2*pi */
  real(dp), parameter :: i2_32m1 = 2.328306437080797e-10 !/* = 1/(2**32 - 1) */
  real(dp), parameter :: d2_32 = 4294967296.d0 ! 2**32
  real(dp), parameter :: i2_30 = 9.31322574615479e-10  ! 2**(-30)

  !---------------------------------------------------------------
  ! Jason Blevins.  Default seed vector
  !---------------------------------------------------------------
  integer, parameter ::  default_seedJB(4) &
       = (/521288629, 362436069, 16163801, 1131199299/)

  !---------------------------------------------------------------
  ! Wichmann-Hill. Default x,y,z
  !---------------------------------------------------------------
  integer, parameter :: default_xyz(3) = (/30269, 30307, 30323/)

  !---------------------------------------------------------------
  ! Mersenne Twister.  Period parameters
  !---------------------------------------------------------------
  integer, parameter :: n_mt = 624
  integer, parameter :: n1_mt = n_mt+1
  integer, parameter :: mata = -1727483681

  !---------------------------------------------------------------
  ! Marsaglia-MultiCarry.  Default seeds x,y,z,w.
  !---------------------------------------------------------------
  integer, parameter :: default_xyzw32(4) &
       = (/123456789, 362436069, 521288629, 916191069/)
  integer(kind=8), parameter :: default_xyzw64(4) &
       = (/1234567890987654321_8, 362436362436362436_8, &
       1066149217761810_8, 123456123456123456_8/)

  !---------------------------------------------------------------
  ! Knuth
  !---------------------------------------------------------------
  integer, parameter :: mm_kn = 1073741824 !2**30
  integer, parameter :: kk_kn = 100  ! the long lag
  integer, parameter :: ll_kn = 37   ! the short lag

  !---------------------------------------------------------------
  !  L'Ecuyer's 1999 random number generator. 64-bits
  !---------------------------------------------------------------
  integer(kind=8), parameter, dimension(5) :: default_seedLe64 = &
       (/153587801,  -759022222,  -759022222, -1718083407, -123456789/)

  type :: rng_t
     integer :: type = 2
     logical :: initialize = .true.

     !---------------------------------------
     ! Jason Blevins. Wichmann-Hill. kiss32
     !---------------------------------------
     integer, dimension(4) :: state

     !----------------------------------------
     !  Mersenne Twister: mt, mti and mag01
     !----------------------------------------
     integer, allocatable :: mt(:)
     integer :: mti = n1_mt   ! mti == N+1 means mt[N] is not initialized
     integer, dimension(0:1) :: mag01 = (/0,mata/)

     !---------------------------------------
     ! Marsaglia-MultiCarry (kiss64). 
     !---------------------------------------
     integer(KIND=8), dimension(5) :: state64

     !----------------------------------------
     !  Knuth
     !----------------------------------------
     integer, allocatable :: randx(:)
     integer :: Kt_pos

     !----------------------------------------
     !   random_standard_gamma
     !----------------------------------------
     real (dp) :: aa = 0, aaa = 0, b, c, d, s, s2, si, q0
  end type rng_t

  interface rexp
     module procedure random_standard_exponential
  end interface rexp

  interface rnorm
     module procedure random_standard_norm
  end interface rnorm


contains

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
  ! 
  !  Subroutines to generate uniform random numbers.
  !  Each subroutine contains a description and references to the original code.
  !  Following the same idea as in rng_seed and rng_uniform by  Jason Blevins, 
  !  all other subroutines were modified to accept as input/output the  
  !  type(rng_t) variable so they can be used in parallel computing.
  !  
  ! Taiane  S. Prass
  ! Porto Alegre
  ! February, 2014
  ! Last revision: March, 2020
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !  http://jblevins.org/log/openmp 
  !  Jason Blevins
  !  Fayetteville, January 27, 2009
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine rng_seed_Blevins(self, seed)
    implicit none
    type(rng_t), intent(inout) :: self
    integer, intent(in) :: seed(4)
    self%state(1:4) = seed
  end subroutine rng_seed_Blevins

  function rng_uniform_Blevins(self) result(u)
    implicit none
    type(rng_t), intent(inout) :: self
    real (dp) :: u
    integer :: imz

    if(self%initialize) then
       call rng_seed_Blevins(self,default_seedJB)
       self%initialize = .false.
    end if

    imz = self%state(1) - self%state(3)
    if (imz < 0) imz = imz + 2147483579
    self%state(1) = self%state(2)
    self%state(2) = self%state(3)
    self%state(3) = imz
    self%state(4) = 69069 * self%state(4) + 1013904243
    imz = imz + self%state(4)
    u = 0.5d0 + 0.23283064d-9 * imz
    return
  end function rng_uniform_Blevins


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !  Wichmann-Hill algorithm  (same as WICHMANN_HILL from R)
  !
  !  Reference:
  !
  !    Brian Wichman, David Hill,
  !    Algorithm AS 183: An Efficient and Portable Pseudo-Random
  !    Number Generator,
  !    Applied Statistics,
  !    Volume 31, Number 2, 1982, pages 188-190.
  !
  !  Discussion:
  !
  !    This function returns a pseudo-random number rectangularly distributed
  !    between 0 and 1.   The cycle length is 6.95E+12.  (See page 123
  !    of Applied Statistics (1984) volume 33), not as claimed in the
  !    original article.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    08 July 2008
  !
  !  Author:
  !
  !    Original FORTRAN77 original version by Brian Wichman, David Hill.
  !    FORTRAN90 version by John Burkardt.
  !  Parameters:
  !
  !    Input/output, integer ( kind = 4 ) S1, S2, S3, three values used as the
  !    seed for the sequence.  These values should be positive
  !    integers between 1 and 30,000.
  !
  !    Output, real ( kind = 4 ) R4_RANDOM, the next value in the sequence.
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Further modified by Taiane Schaedler Prass, 
  ! now accepts the arguments self and seed
  ! For compatibility with other subroutines uses "dp" instead of "kind = 4" 
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine rng_seed_wh(self, seed)
    implicit none
    type(rng_t), intent(inout) :: self
    integer, intent(in)   :: seed(3)
    self%state(1:3) = seed
    return
  end subroutine rng_seed_wh

  function rng_uniform_wh(self) result(fn_val)
    implicit none
    type(rng_t), intent(inout) :: self
    real(dp) :: fn_val

    if(self%initialize) call rng_seed_wh(self, default_xyz)

    self%state(1) = mod(171*self%state(1), 30269)
    self%state(2) = mod(172*self%state(2), 30307)
    self%state(3) = mod(170*self%state(3), 30323)
    fn_val = mod(dble(self%state(1))/30269.0d0 + &
         dble(self%state(2))/30307.0d0 + &
         dble(self%state(3))/30323.0d0, 1.d0)

    return
  end function rng_uniform_wh


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !   Mersenne Twister algorithm.  (same as MERSENNE_TWISTER from R)
  !
  !   https://jblevins.org/mirror/amiller/mt19937.f90   
  !   https://jblevins.org/mirror/amiller/mt19937a.f90
  !   Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
  !   Code converted using TO_F90 by Alan Miller.   Date: 1999-11-26  Time: 17:09:23
  !   Latest revision - 5 February 2002
  !
  !   Modified to accept self and seed arguments.
  !   by Taiane S. Prass
  !   Porto Alegre, 2012
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine rng_seed_sgrnd(self, seed)
    implicit none
    type(rng_t), intent(inout) :: self
    integer, intent(in)   :: seed
    real(dp), parameter :: two31 = 2147483648.d0 ! 2**31
    real(dp), parameter :: i2_31 = 4.65661287307739e-10  ! 1/2**31
    real (dp)  :: temp
    integer    :: itemp, itemp2, mti

    if(allocated(self%mt)) deallocate(self%mt)
    allocate(self%mt(0:(n_mt-1)))
    self%mt = 0

    !----------------------------------------------------------------------------------
    !  setting initial seeds to mt[N] using the generator Line 25 of Table 1 in
    !  [KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102]
    !----------------------------------------------------------------------------------
    self%mt(0)= iand(seed, -1)
    do  mti = 1, (n_mt-1)
       !-------------------------------------------------------------
       ! The following code in this loop is equivalent to:
       ! the single line of code:
       !            mt(mti) = iand(69069 * mt(mti-1), -1)
       ! The code here is used instead to prevent integer overflow.
       !-------------------------------------------------------------
       temp = 69069.d0*dble(self%mt(mti-1))
       itemp = int(mod(temp,two31))
       itemp2 = int(temp*i2_31)
       if (mod(itemp2,2).ne.0) then
          if (itemp.gt.0) then
             itemp = int(itemp - two31)
          else
             itemp = int(itemp + two31)
          endif
       endif
       self%mt(mti) = itemp       
    end do
    self%mti = mti
    return
  end subroutine rng_seed_sgrnd

  function rng_uniform_Mersenne(self) result(fn_val)
    implicit none
    type(rng_t), intent(inout) :: self
    real (dp) :: fn_val

    integer, parameter :: lmask =  2147483647        !  least significant r bits
    !integer, parameter :: umask = -lmask - 1  !  most significant w-r bits
    integer, parameter :: tmaskb= -1658038656, tmaskc= -272236544
    integer, parameter :: m_mt = 397, seed = 4357
    integer  :: k, y, umask

    umask = ISHFT(1073741824,1)
    !generate N words at one time
    if(self%initialize  .or. self%mti == n_mt+1) then
       !  if sgrnd() has not been called,
       !  a default initial seed is used
       call rng_seed_sgrnd(self, seed) 
       self%initialize = .false.
    end if

    if(self%mti >= n_mt) then      
       do k = 0, (n_mt - m_mt - 1)
          y = ior(iand(self%mt(k),umask), iand(self%mt(k+1),lmask))
          self%mt(k) = ieor(ieor(self%mt(k+m_mt), ishft(y,-1)),self%mag01(iand(y,1)))
       end do
       do  k = (n_mt - m_mt), (n_mt - 2)
          y = ior(iand(self%mt(k),umask), iand(self%mt(k+1),lmask))
          self%mt(k) = ieor(ieor(self%mt(k+(m_mt-n_mt)), ishft(y,-1)),self%mag01(iand(y,1)))
       end do
       y = ior(iand(self%mt(n_mt-1),umask), iand(self%mt(0),lmask))
       self%mt(n_mt-1) = ieor(ieor(self%mt(m_mt-1), ishft(y,-1)),self%mag01(iand(y,1)))
       self%mti = 0
    end if

    y = self%mt(self%mti)
    self%mti = self%mti + 1
    y = ieor(y, ishft(y,-11))
    y = ieor(y, iand(ishft(y,7),tmaskb))
    y = ieor(y, iand(ishft(y,15),tmaskc))
    y = ieor(y, ishft(y,-18))

    if(y < 0) then
       fn_val = (dble(y) + d2_32)*i2_32m1 !1 - y/2**32
    else
       fn_val = dble(y)*i2_32m1  ! y/2**32
    end if
    return
  end function rng_uniform_Mersenne


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! public domain code
  ! made available from
  ! http://www.fortran.com/kiss.f90
  !
  ! The  KISS (Keep It Simple Stupid) random number generator. Combines:
  ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
  ! (2) A 3-shift shift-register generator, period 2^32-1,
  ! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
  !  Overall period>2^123;  Default seeds x,y,z,w.
  !  Set your own seeds with statement i=kisset(ix,iy,iz,iw).            !
  !
  !
  !  Modified to accept self and seed arguments.
  !  by Taiane S. Prass
  !   Porto Alegre, 2012
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine rng_seed_kiss32(self, seed)
    implicit none
    type(rng_t), intent(inout) :: self
    integer, intent(in) :: seed(4)
    self%state(1:4) = seed
    return
  end subroutine rng_seed_kiss32

  function m_kiss32(k, n) result(m)
    implicit none
    integer, intent(in) :: k,n
    integer :: m
    m = ieor (k, ishft (k, n) )
    return
  end function m_kiss32

  function  rng_uniform_kiss32(self) result(unif_rand)
    implicit none
    type(rng_t), intent(inout) :: self
    real (dp) :: unif_rand
    integer :: kiss
    real(dp), parameter :: i2_32m1 = 2.328306437080797e-10 !/* = 1/(2**32 - 1) */
    real(dp), parameter :: d2_32 = 4294967296.d0 ! 2**32

    if(self%initialize) then
       call rng_seed_kiss32(self, default_xyzw32) 
       self%initialize = .false.
    end if

    self%state(1) = 69069 * self%state(1) + 1327217885
    self%state(2) = m_kiss32 (m_kiss32 (m_kiss32 (self%state(2), 13), - 17), 5)
    self%state(3) = 18000 * iand (self%state(3), 65535) + ishft (self%state(3), - 16)
    self%state(4) = 30903 * iand (self%state(4), 65535) + ishft (self%state(4), - 16)
    kiss = self%state(1) + self%state(2) + ishft (self%state(3), 16) + self%state(4)

    if(kiss < 0) then
       unif_rand = (kiss + d2_32)*i2_32m1
    else
       unif_rand = kiss*i2_32m1
    end if

    return
  end function rng_uniform_kiss32

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! http://lgge.osug.fr/meom/pages-perso/brankart/Outils/mod_kiss.f90
  ! by Jean-Michel Brankart
  !
  ! Based on the
  ! 64-bit KISS (Keep It Simple Stupid) random number generator
  ! distributed by George Marsaglia :
  ! http://groups.google.com/group/comp.lang.fortran/
  !        browse_thread/thread/a85bf5f2a97f5a55
  !
  ! The 64-bit KISS (Keep It Simple Stupid) random number generator.
  ! Components:
  !  (1) Xorshift (XSH), period 2^64-1,
  !  (2) Multiply-with-carry (MWC), period (2^121+2^63-1)
  !  (3) Congruential generator (CNG), period 2^64.
  ! Overall period:
  !  (2^250+2^192+2^64-2^186-2^129)/6 ~= 2^(247.42) or 10^(74.48)
  ! Set your own seeds with statement "CALL kiss_seed(ix,iy,iz,iw)".
  !
  !  Modified to accept self and seed arguments.
  !  by Taiane S. Prass
  !   Porto Alegre, 2012
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine rng_seed_kiss64(self, seed)
    implicit none
    type(rng_t), intent(inout) :: self
    integer, intent(in) :: seed(4)
    integer :: i
    do i = 1,4
       self%state64(i) = default_xyzw64(i)+int(seed(i),8)
    end do
    return
  end subroutine rng_seed_kiss64

  function m_kiss64(k, n) result(m)
    implicit none
    integer(kind=8), intent(in) :: k,n
    integer(kind=8) :: m
    m = ieor (k, ishft (k, n) )
    return
  end function m_kiss64

  function rng_uniform_kiss64(self) result(fn_val)
    implicit none
    type(rng_t), intent(inout) :: self
    integer(kind=8) ::  kiss, t
    real (dp) :: fn_val
    real (dp), parameter :: huge64=9223372036854775808.0d0  ! 2**63

    if(self%initialize) then
       call rng_seed_kiss64(self, (/0,0,0,0/)) 
       self%initialize = .false.
    end if

    t = ishft(self%state64(1), 58) + self%state64(4)
    if (ishft(self%state64(1), -63)  == ishft(t, -63)) then
       self%state64(4) = ishft(self%state64(1),-6) + ishft(self%state64(1), -63)
    else
       self%state64(4) = ishft(self%state64(1),-6) + 1 - ishft(self%state64(1)+t, -63)
    end if
    self%state64(1) = t + self%state64(1)
    self%state64(2) = m_kiss64( m_kiss64( m_kiss64(self%state64(2),13_8), -17_8 ), 43_8 )
    self%state64(3) = 6906969069_8 * self%state64(3) + 1234567_8
    kiss = self%state64(1) + self%state64(2) + self%state64(3)

    !-----------------------------------------------------------
    ! Real random numbers with uniform distribution in [0,1]
    !-----------------------------------------------------------
    fn_val = 0.5d0 * ( 1.d0 + real(kiss, 8)/huge64 )
    return
  end function rng_uniform_kiss64

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! https://jblevins.org/mirror/amiller/rand3.f90
  !
  !   Random Number Generator based on Knuth's program 
  !     https://www-cs-faculty.stanford.edu/~knuth/programs/frng.f
  !
  !  Here:
  !   rng_array replaces the f77 version ran_array
  !   rng_seed_rnstrt replaces rnstrt
  !   rng_uniform_knuth combines the subroutines to create uniform random numbers
  !
  !   Reference: Seminumerical Algorithms, 3rd edition, Section 3.6
  !   
  !   Note: includes the modifications made in the 9th printing (2002),
  !
  ! Author: Steve Kifowit
  ! http://ourworld.compuserve.com/homepages/steve_kifowit
  ! with modifications by Alan Miller to rnarry and rnstrt based upon
  ! Knuth's code.
  !
  !  Modified to accept self and seed arguments.
  !  by Taiane S. Prass
  !   Porto Alegre, 2012.
  !   Last revision, 2019
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine rng_array(aa, n, self)
    !--------------------------------------------------
    ! creates an array using the recurrency
    !  X(j) = [X(j-100) - X(j-37)] mod 2**30
    !--------------------------------------------------
    implicit none
    type(rng_t), intent(inout) :: self
    integer, intent(in) :: n
    integer :: aa(n)   ! auxiliar array
    integer :: j

    !------------------------------------------------
    !  aa(j) = X(j), j = 1,...,n   n = quality level
    !  aa(j) = aa(j-100) - aa(j-37), j = 101,...,n           
    !------------------------------------------------
    do  j = 1,kk_kn
       aa(j) = self%randx(j)
    end do
    do  j = kk_kn+1,n
       aa(j) = aa(j-kk_kn) - aa(j-ll_kn)
       if(aa(j) < 0) aa(j) = aa(j) +  mm_kn
    end do

    !------------------------------------------------------
    !  randx(j) = aa(n+j-100) - aa(n+j-37), j = 1,...,37
    !  randx(j) = aa(n+j-100) - X(j-37), j = 38,...,100
    !------------------------------------------------------
    do  j = 1,ll_kn
       self%randx(j) = aa(n+j-kk_kn) - aa(n+j-ll_kn)
       if(self%randx(j)  < 0) self%randx(j) = self%randx(j) +  mm_kn
    end do
    do  j = ll_kn+1, kk_kn
       self%randx(j) = aa(n+j-kk_kn) - self%randx(j-ll_kn)
       if(self%randx(j) < 0) self%randx(j) = self%randx(j) +  mm_kn
    end do

    return
  end subroutine rng_array

  subroutine rng_seed_rnstrt(self, seed)
    implicit none
    type(rng_t), intent(inout) :: self
    integer, intent(in) :: seed
    integer, parameter :: tt = 70  ! guaranteed separation between streams 
    integer :: x(2*kk_kn-1) ! the preparation buffer 
    integer :: j, ss, t, i, sseed, kmax

    kmax = 2*kk_kn-1    
    !------------------------------------------------------------
    ! allocating the variable to save the values for future call
    !------------------------------------------------------------
    if(allocated(self%randx)) deallocate(self%randx)
    allocate(self%randx(kk_kn))    
    self%randx = 0

    !-----------------------------------------------
    ! Make the seed even
    !-----------------------------------------------     
    sseed = mod(seed, mm_kn)
    ss = sseed + 2 - mod(sseed,2)

    do  j = 1,kk_kn
       x(j) = ss
       ss = ss + ss
       if (ss >= mm_kn) ss = ss - mm_kn + 2
    end do
    ! make (only) x(2) odd
    x(2) = x(2)+1

    ss = sseed
    t = tt-1
10  continue
    !-----------------------------------------------
    ! Loop to fill x
    !-----------------------------------------------
    do j = kk_kn,2,-1
       x(j+j-1)=x(j)
       x(j+j-2)=0
    end do

    do j = kmax,(kk_kn+1),-1     
       x(j-(kk_kn-ll_kn)) =x(j-(kk_kn-ll_kn))-x(j)
       if(x(j-(kk_kn-ll_kn)) < 0)  x(j-(kk_kn-ll_kn)) =  x(j-(kk_kn-ll_kn)) + mm_kn
       x(j-kk_kn) = x(j-kk_kn)-x(j)
       if(x(j-kk_kn) < 0) x(j-kk_kn) = x(j-kk_kn) + mm_kn
    end do

    if (mod(ss,2) == 1) then
       do  j = kk_kn,1,-1
          x(j+1)=x(j)
       end do
       x(1) = x(kk_kn+1)
       x(ll_kn+1) = x(ll_kn+1)-x(kk_kn+1)
       if(x(ll_kn+1) < 0 )  x(ll_kn+1) =  x(ll_kn+1) + mm_kn
    end if

    if (ss /= 0) then
       ss=ss/2
    else
       t = t-1
    end if
    if (t > 0) go to 10

    !-----------------------------------------------
    ! saving initial values
    !-----------------------------------------------
    do  j = 1,ll_kn         ! randx(64:100)
       self%randx(j+kk_kn-ll_kn) = x(j)
    end do
    do j = (ll_kn+1),kk_kn  ! randx(1:63)       
       self%randx(j-ll_kn) = x(j)
    end do

    !-----------------------------------------------
    ! warm things up
    !-----------------------------------------------
    do i = 1,10
       call rng_array(x, kmax, self)
    end do

    self%Kt_pos = 1
    return
  end subroutine rng_seed_rnstrt


  function rng_uniform_knuth(self) result(fn_val)
    !-------------------------------------------
    !  Knuth's random number generator
    !-------------------------------------------
    implicit none
    type(rng_t), intent(inout) :: self
    real (dp) :: fn_val
    integer, parameter :: ql_kn = 1009  ! quality level
    integer, parameter :: seed = 314159
    integer :: aa(ql_kn)

    if(self%initialize) then
       call rng_seed_rnstrt(self,seed)
       self%initialize = .false.
    end if

    ! If all 100 values where used, restart the cycle 
    if(self%Kt_pos > 100) then
       ! rng_array_cycle using quality level = 1009
       call rng_array(aa, ql_kn, self)       
       self%Kt_pos = 1
    end if
    fn_val = scale(dble(self%randx(self%Kt_pos)), -30)
    self%Kt_pos = self%Kt_pos  + 1     

    return
  end function rng_uniform_knuth

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! L'Ecuyer's 1999 random number generator.
  ! Fortran version by Alan.Miller @ vic.cmis.csiro.au
  ! This version is for 64-bit integers and assumes that KIND=8 identifies them.
  ! Latest revision - 12 January 2001
  ! https://jblevins.org/mirror/amiller/lfsr258.f90
  !
  ! Generates a random number between 0 and 1.  Translated from C function in:
  ! Reference:
  ! L'Ecuyer, P. (1999) `Tables of maximally equidistributed combined LFSR
  ! generators', Math. of Comput., 68, 261-269.
  ! The cycle length is claimed to be about 2^(258) or about 4.6 x 10^77.
  ! Actually - (2^63 - 1).(2^55 - 1).(2^52 - 1).(2^47 - 1).(2^41 - 1)
  !
  !  Modified to accept self and seed arguments.
  !  by Taiane S. Prass
  !   Porto Alegre, 2012
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine rng_seed_lfsr258(self, seed)
    implicit none
    type(rng_t), intent(inout) :: self
    integer, intent(in) :: seed(5)
    integer :: i
    do i = 1,5       
       self%state64(1:5) = default_seedLe64 + int(seed(i), 8)
    end do
    if (iand(self%state64(1), -int(2,8)) == 0) self%state64(1) = self%state64(1) - 8388607
    if (iand(self%state64(2), -int(512,8)) == 0) self%state64(2) = self%state64(2) - 8388607
    if (iand(self%state64(3), -int(4096,8)) == 0) self%state64(3) = self%state64(3) - 8388607
    if (iand(self%state64(4), -int(131072,8)) == 0) self%state64(4) = self%state64(4) - 8388607
    if (iand(self%state64(5),-int(8388608,8)) == 0) self%state64(5) = self%state64(5) - 8388607
    return
  end subroutine rng_seed_lfsr258

  function rng_uniform_Le(self) result(random_numb)
    implicit none
    type(rng_t), intent(inout) :: self
    real (dp) :: random_numb
    integer (kind = 8)  :: b, dummy

    if(self%initialize) then
       call rng_seed_lfsr258(self,(/0,0,0,0,0/))
       self%initialize = .false.
    end if

    !----------------------------------------------------------------
    ! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
    !      to the left if j > 0, otherwise to the right.
    !----------------------------------------------------------------
    b  = ishft( ieor( ishft(self%state64(1),1), self%state64(1)), -53)
    self%state64(1) = ieor( ishft( iand(self%state64(1),-int(2,8)), 10), b)
    b  = ishft( ieor( ishft(self%state64(2),24), self%state64(2)), -50)
    self%state64(2) = ieor( ishft( iand(self%state64(2),-int(512,8)), 5), b)
    b  = ishft( ieor( ishft(self%state64(3),3), self%state64(3)), -23)
    self%state64(3) = ieor( ishft( iand(self%state64(3),-int(4096,8)), 29), b)
    b  = ishft( ieor( ishft(self%state64(4),5), self%state64(4)), -24)
    self%state64(4) = ieor( ishft( iand(self%state64(4),-int(131072,8)), 23), b)
    b  = ishft( ieor( ishft(self%state64(5),3), self%state64(5)), -33)
    self%state64(5) = ieor( ishft( iand(self%state64(5),-int(8388608,8)), 8), b)
    !----------------------------------------------------------------
    ! The constant below is the reciprocal of (2^64 - 1)
    !----------------------------------------------------------------
    dummy = ieor( ieor( ieor( ieor(self%state64(1),self%state64(2)), &
         self%state64(3)), self%state64(4)), self%state64(5))
    random_numb =  dble(dummy) * 5.4210108624275221e-20 + 0.5d0
    return
  end function rng_uniform_Le

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !   Main surboutines: rng_seed and rng_uniform
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
  subroutine rng_seed(self, ns, seed, type)
    implicit none
    type(rng_t), intent(inout) :: self
    integer, intent(in) :: ns
    integer, intent(in) :: seed(ns), type
    integer, allocatable :: dummy(:)

    self%type = type
    self%initialize = .false.

    select case(self%type)
    case (0)
       !----------------------------
       !   Jason Blevins. ns = 4
       !----------------------------
       allocate(dummy(4))
       if(ns < 4) then       
          dummy(1:ns) = seed
          dummy((ns+1):4) = default_seedJB((ns+1):4)
       else
          dummy = seed(1:4)
       end if
       call rng_seed_Blevins(self,dummy)
    case (1)
       !-------------------------
       !  Wichmann-Hill. ns = 3
       !-------------------------
       allocate(dummy(3))
       if(ns < 3) then       
          dummy(1:ns) = seed
          dummy((ns+1):3) = default_xyz((ns+1):3)
       else
          dummy = seed(1:3)
       end if
       call rng_seed_wh(self,dummy)      
    case (2)
       !-----------------------------
       !   Mersenne Twister. ns = 1
       !-----------------------------
       call  rng_seed_sgrnd(self,seed(1))
    case (3)
       !----------------------------------------
       ! Marsaglia-MultiCarry (kiss 32). ns = 4
       !----------------------------------------
       allocate(dummy(4))
       if(ns < 4) then       
          dummy(1:ns) = seed
          dummy((ns+1):4) = default_xyzw32((ns+1):4)
       else
          dummy = seed(1:4)
       end if
       call rng_seed_kiss32(self,dummy)
    case (4)
       !----------------------------------------
       ! Marsaglia-MultiCarry (kiss 64) ns = 4
       !----------------------------------------
       allocate(dummy(4))
       if(ns < 4) then       
          dummy(1:ns) = seed
          dummy((ns+1):4) = 0
       else
          dummy = seed(1:4)
       end if
       call rng_seed_kiss64(self,dummy)
    case (5)
       !---------------------------------
       ! Knuth (2002). ns = 1
       !---------------------------------
       call rng_seed_rnstrt(self,seed(1))
    case (6)
       !-----------------------------------
       ! L'Ecuyer's (1999, 64-bits). ns = 5
       !-----------------------------------
       allocate(dummy(5))
       if(ns < 5) then       
          dummy(1:ns) = seed
          dummy((ns+1):5) = 0
       else
          dummy = seed(1:5)
       end if
       call rng_seed_lfsr258(self,dummy)
    end select
    return
  end subroutine rng_seed


  function rng_uniform(self) result(fn_val)
    implicit none
    type(rng_t), intent(inout) :: self
    real (dp) :: fn_val

    select case (self%type)
    case (0) ! Jason Blevins
       fn_val = rng_uniform_Blevins(self)
    case (1) ! Wichmann-Hill
       fn_val = rng_uniform_wh(self)
    case (2) ! Mersenne Twister
       fn_val = rng_uniform_mersenne(self)
    case (3) ! Marsaglia-MultiCarry (kiss 32)
       fn_val = rng_uniform_kiss32(self)
    case (4) ! Marsaglia-MultiCarry (kiss 64)
       fn_val = rng_uniform_kiss64(self)
    case (5) ! Knuth (2002)
       fn_val = rng_uniform_Knuth(self)
    case (6) ! L'Ecuyer's 1999 (64-bits)
       fn_val = rng_uniform_Le(self)
    case default
       fn_val = rng_uniform_mersenne(self)
    end select

    return
  end function rng_uniform


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
  !
  !    Beta distribution
  !  
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
  function rbeta(npar, par, rng) result(fn)
    !
    ! Adapted for compatibility with other distributions
    !
    implicit none
    type(rng_t), intent(inout) :: rng
    integer, intent(in) :: npar   ! must be 2
    real(dp), intent(in) :: par(npar)      
    real(dp) :: fn
    fn = rbeta_default(par(1)*par(2), (1.d0 - par(1))*par(2), rng)
    return
  end function rbeta


  function rbeta_default(a, b, rng) result(fn)
    !
    !  Inspired in the function in R. 
    !
    implicit none
    real(dp), intent(in) :: a, b
    type(rng_t), intent(inout) :: rng
    real(dp) :: fn

    if(a < 0.d0 .or. b < 0.d0) then
       ! Erro
       fn = 999.d0
       return
    end if

    if(a > Huge(1.d0) .and. b > Huge(1.d0)) then
       ! a = b = Inf: all mass at 1/2
       fn = 0.5d0
       return
    end if

    if(a == 0.0d0 .and. b == 0.0d0) then
       ! a = b = 0: point mass 1/2 at each of {0,1}
       fn = 1.d0
       if(rng_uniform(rng) < 0.5d0) fn = 0.d0
       return
    end if

    if(a == 0.0d0) then
       fn = 0.d0
       return
    end if

    if(b == 0.0d0) then
       fn = 1.d0
       return
    end if

    fn = random_beta(a,b,rng)
    return
  end function rbeta_default

  FUNCTION random_beta(aa, bb, rng) RESULT(fn_val)
    !---------------------------------------------------------
    !
    ! function to generate a beta random variable
    ! adapted to be used in parallel
    !
    !---------------------------------------------------------
    ! https://jblevins.org/mirror/amiller/random.f90
    !
    ! Adapted from Fortran 77 code from the book:
    !     Dagpunar, J. 'Principles of random variate generation'
    !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
    !
    ! Author: Alan Miller
    ! e-mail: amiller @ bigpond.net.au
    !
    ! FUNCTION GENERATES A RANDOM VARIATE IN [0,1]
    ! FROM A BETA DISTRIBUTION WITH DENSITY
    ! PROPORTIONAL TO BETA**(AA-1) * (1-BETA)**(BB-1).
    ! USING CHENG'S LOG LOGISTIC METHOD.
    !
    !     AA = SHAPE PARAMETER FROM DISTRIBUTION (0 < REAL)
    !     BB = SHAPE PARAMETER FROM DISTRIBUTION (0 < REAL)
    !
    REAL(dp), INTENT(IN) :: aa, bb
    type(rng_t), intent(inout) :: rng
    REAL(dp) :: fn_val

    !     Local variables
    REAL(dp), PARAMETER  :: aln4 = 1.3862943611198906, vlarge = HUGE(1.d0), vsmall = TINY(1.d0)
    real(dp), parameter :: one = 1.d0, zero = 0.d0, two = 2.d0
    REAL(dp) :: a, b, g, r, s, x, y, z
    REAL(dp) :: d, f, h, t, c
    LOGICAL :: swap

    IF (aa <= zero .OR. bb <= zero) THEN
       !WRITE(*, *) 'IMPERMISSIBLE SHAPE PARAMETER VALUE(S)'
       !STOP
       fn_val = 999.d0
       return
    END IF

    a = aa
    b = bb
    swap = b > a
    IF (swap) THEN
       g = b
       b = a
       a = g
    END IF
    d = a/b
    f = a+b
    IF (b > one) THEN
       h = SQRT((two*a*b - f)/(f - two))
       t = one
    ELSE
       h = b
       t = one/(one + (a/(vlarge*b))**b)
    END IF
    c = a+h

    fn_val = zero

    DO
       r = rng_uniform(rng)
       x = rng_uniform(rng)
       s = r*r*x
       IF (r < vsmall .OR. s <= zero) CYCLE
       IF (r < t) THEN
          x = LOG(r/(one - r))/h
          y = d*EXP(x)
          z = c*x + f*LOG((one + d)/(one + y)) - aln4
          IF (s - one > z) THEN
             IF (s - s*z > one) CYCLE
             IF (LOG(s) > z) CYCLE
          END IF
          fn_val = y/(one + y)
       ELSE
          IF (4.d0*s > (one + one/d)**f) CYCLE
          fn_val = one
       END IF
       EXIT
    END DO

    IF (swap) fn_val = one - fn_val
    RETURN
  END FUNCTION random_beta


  function dbeta_default(x,a,b,give_log) result(fn_val)
    !-----------------------------------------------------------------------
    !               evaluation of x**(a-1)*(1-x)**(b-1)/beta(a,b)
    !-----------------------------------------------------------------------
    implicit none
    real (dp),  intent(in)   :: x, a, b
    logical, intent(in) :: give_log
    real(dp) :: fn_val, y

    fn_val = 0.0d0
    y =  1.0d0 - x
    if (x <= 0.0d0 .or. y <= 0.0d0) return
    fn_val = log(brcomp(a,b,x,y))- log(x) -log(y)
    if(.not. give_log) fn_val = exp(fn_val)      
    return
  end function dbeta_default

  function dbeta(x, npar, par, give_log) result(fn_val)
    implicit none
    integer, intent(in) :: npar
    real (dp),  intent(in)   :: x, par(npar)
    logical, intent(in) :: give_log
    real(dp) :: fn_val 
    fn_val = dbeta_default(x,par(1)*par(2),(1.d0 - par(1))*par(2), give_log)      
    return
  end function dbeta


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
  !
  !     Kumaraswamy distribution
  !  
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function rkuma_default(delta, nu, a, b, rng) result(fn)
    !----------------------------------------------------
    !
    !  inversion method for randon number generation
    !
    ! nu/varphi > 0 
    ! delta > 0
    ! a  = lower limit
    ! b  = upper limit
    !----------------------------------------------------
    implicit none
    real(dp), intent(in) :: delta, nu,  a, b
    type(rng_t), intent(inout) :: rng
    real(dp) :: u, fn

    u = rng_uniform(rng)
    ! To avoid overflow/underflow
    u = (1.d0/delta)*log(1.d0 - u)
    u = exp(u)                       !z = (1-u)**(1/delta)
    u = 1.d0/nu*log(1.d0 - u)   
    fn = exp(u)                      !x = (1 - z)**1/nu
    fn = a + (b-a)*fn                !y = a + (b-a)*x

    return
  end function rkuma_default

  function rkuma(npar, par, rng) result(fn)
    ! par = (mu, nu, rho, a, b)
    implicit none
    integer, intent(in) :: npar   ! must be 5
    real(dp), intent(in) :: par(npar)
    type(rng_t), intent(inout) :: rng
    real(dp) :: fn, delta, mu01

    mu01 = (par(1) - par(4)) / (par(5) - par(4))
    delta = log(1.d0 - par(3)) / log( 1.d0 - mu01**par(2))
    fn = rkuma_default(delta, par(2), par(4), par(5), rng)
    return
  end function rkuma


  function dkuma(y, npar, par, give_log) result(fn)
    !-----------------------------------------------------
    ! Density funtion - Kumaraswamy distribution
    ! a, b real numbers 
    !
    ! mu = rho-quantile, 0 < rho < 1
    ! 
    ! par = c(mu, nu, rho, a, b)
    !
    ! Implemented by Taiane S. Prass
    ! April, 2020
    !
    ! Last revision: April, 2021 
    !    Replace original dkuma by dkuma and dkuma_default
    !-----------------------------------------------------
    implicit none
    integer, intent(in):: npar
    real(dp), intent(in) :: y, par(npar)
    logical, intent(in) :: give_log      
    real(dp) :: fn, mu01, dt

    mu01 = (par(1) - par(4)) / (par(5) - par(4))
    dt = log(1.d0 - par(3)) / log(1.d0 - mu01**par(2))       
    fn = dkuma_default(y, dt, par(2), par(4), par(5), give_log)
    return
  end function dkuma


  function dkuma_default(y, delta, nu, a, b, give_log) result(fn)
    implicit none
    real(dp), intent(in) :: y, delta, nu, a, b
    logical, intent(in) :: give_log      
    real(dp) :: fn

    fn = log(nu) + log(delta) - log(b-a) 
    fn = fn + (nu - 1.d0)*(log(y-a) - log(b-a))
    fn = fn + (delta - 1.d0)*log(1.d0 - ((y - a) / (b - a))**nu)   
    if(.not. give_log) fn = exp(fn)
    return
  end function dkuma_default


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
  !
  !     Exponential distribution
  !  
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function random_standard_exponential(rng)
    implicit none
    type(rng_t), intent(inout) :: rng
    real (dp) :: random_standard_exponential
    real (dp) :: a, u, umin, ustar
    integer :: i
    real (dp), parameter :: q(8) =  (/0.6931472, 0.9333737, 0.9888778, &
         0.9984959, 0.9998293,  0.9999833, 0.9999986, 0.9999999/)
    a = 0.0
    u = rng_uniform(rng)
    go to 20
10  continue
    a = a + q(1)
20  continue
    u = u + u
    if (u < 1.0) go to 10
    u = u - 1.0
    if (u <= q(1)) then
       random_standard_exponential = a + u
       return
    end if
    i = 1
    ustar = rng_uniform(rng)
    umin = ustar
    ustar = rng_uniform(rng)
    umin = min(ustar,umin)
    i = i + 1
    do while (u > q(i))
       ustar = rng_uniform(rng)
       umin = min(ustar,umin)
       i = i + 1
    end do
    random_standard_exponential = a + umin*q(1)
    return
  end function random_standard_exponential


  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
  !
  !     Normal distribution
  !  
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function random_standard_norm(rng) result(fn)
    !---------------------------------------------------------
    !  Same as rnorm from R
    !---------------------------------------------------------
    implicit none
    type(rng_t), intent(inout) :: rng
    integer, parameter :: BIG = 134217728 !/* 2^27 */
    !/* unif_rand() alone is not of high enough precision */
    real(dp) :: u1, fn

    u1 = rng_uniform(rng)
    u1 = BIG*u1 + rng_uniform(rng)
    fn = qnorm(u1/BIG, 0._dp, 1._dp)
    return 
  end function random_standard_norm

  function qnorm(p, mean, sd) result(q)
    !-------------------------------------
    ! quantile for normal distribution
    !-------------------------------------
    implicit none
    real (dp), intent(in) :: p, mean, sd
    real(dp) :: q
    integer :: ifault
    call standard_qnorm(p,q,ifault)
    q = q*sd + mean
    return
  end function qnorm

  subroutine standard_qnorm (p, normal_dev, ifault)
    !---------------------------------------------------------------
    !
    ! https://jblevins.org/mirror/amiller/as241.f90
    !
    ! ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
    !
    ! Produces the normal deviate Z corresponding to a given lower
    ! tail area of P; Z is accurate to about 1 part in 10**16.
    !
    ! The hash sums below are the sums of the mantissas of the
    ! coefficients.   They are included for use in checking
    ! transcription.
    !
    ! This ELF90-compatible version by Alan Miller - 20 August 1996
    ! N.B. The original algorithm is as a function; this is a subroutine
    !---------------------------------------------------------------
    implicit none
    real (dp), intent(in)   :: p
    integer, intent(out)    :: ifault
    real (dp), intent(out)  :: normal_dev

    ! Coefficients for P close to 0.5
    real (dp), parameter :: a0 = 3.3871328727963666080_dp
    real (dp), parameter :: a1 = 1.3314166789178437745D+2
    real (dp), parameter :: a2 = 1.9715909503065514427D+3
    real (dp), parameter :: a3 = 1.3731693765509461125D+4
    real (dp), parameter :: a4 = 4.5921953931549871457D+4
    real (dp), parameter :: a5 = 6.7265770927008700853D+4
    real (dp), parameter :: a6 = 3.3430575583588128105D+4
    real (dp), parameter :: a7 = 2.5090809287301226727D+3
    real (dp), parameter :: b1 = 4.2313330701600911252D+1
    real (dp), parameter :: b2 = 6.8718700749205790830D+2
    real (dp), parameter :: b3 = 5.3941960214247511077D+3
    real (dp), parameter :: b4 = 2.1213794301586595867D+4
    real (dp), parameter :: b5 = 3.9307895800092710610D+4
    real (dp), parameter :: b6 = 2.8729085735721942674D+4
    real (dp), parameter :: b7 = 5.2264952788528545610D+3
    ! HASH SUM AB    55.8831928806149014439

    ! Coefficients for P not close to 0, 0.5 or 1.
    real (dp), parameter :: c0 = 1.42343711074968357734_dp
    real (dp), parameter :: c1 = 4.63033784615654529590_dp
    real (dp), parameter :: c2 = 5.76949722146069140550_dp
    real (dp), parameter :: c3 = 3.64784832476320460504_dp
    real (dp), parameter :: c4 = 1.27045825245236838258_dp
    real (dp), parameter :: c5 = 2.41780725177450611770D-1
    real (dp), parameter :: c6 = 2.27238449892691845833D-2
    real (dp), parameter :: c7 = 7.74545014278341407640D-4
    real (dp), parameter :: d1 = 2.05319162663775882187_dp
    real (dp), parameter :: d2 = 1.67638483018380384940_dp
    real (dp), parameter :: d3 = 6.89767334985100004550D-1
    real (dp), parameter :: d4 = 1.48103976427480074590D-1
    real (dp), parameter :: d5 = 1.51986665636164571966D-2
    real (dp), parameter :: d6 = 5.47593808499534494600D-4
    real (dp), parameter :: d7 = 1.05075007164441684324D-9
    ! HASH SUM CD    49.33206503301610289036

    ! Coefficients for P near 0 or 1.
    real (dp), parameter :: e0 = 6.65790464350110377720_dp
    real (dp), parameter :: e1 = 5.46378491116411436990_dp
    real (dp), parameter :: e2 = 1.78482653991729133580_dp
    real (dp), parameter  :: e3 = 2.96560571828504891230D-1
    real (dp), parameter :: e4 = 2.65321895265761230930D-2
    real (dp), parameter :: e5 = 1.24266094738807843860D-3
    real (dp), parameter :: e6 = 2.71155556874348757815D-5
    real (dp), parameter :: e7 = 2.01033439929228813265D-7
    real (dp), parameter :: f1 = 5.99832206555887937690D-1
    real (dp), parameter :: f2 = 1.36929880922735805310D-1
    real (dp), parameter :: f3 = 1.48753612908506148525D-2
    real (dp), parameter :: f4 = 7.86869131145613259100D-4
    real (dp), parameter :: f5 = 1.84631831751005468180D-5
    real (dp), parameter :: f6 = 1.42151175831644588870D-7
    real (dp), parameter :: f7 = 2.04426310338993978564D-15
    ! HASH SUM EF    47.52583317549289671629

    real (dp), parameter :: zero = 0._dp, one = 1._dp, half = 0.5_dp, const2 = 1.6_dp
    real (dp), parameter :: split1 = 0.425_dp, split2 = 5._dp, const1 = 0.180625_dp
    real (dp) ::  q, r

    ifault = 0
    q = p - half
    if (abs(q) <= split1) then
       r = const1 - q * q
       normal_dev = q * (((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0) / &
            (((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + one)
       return
    else
       if (q < zero) then
          r = p
       else
          r = one - p
       end if
       if (r <= zero) then
          ifault = 1
          normal_dev = zero
          return
       end if
       r = sqrt(-log(r))
       if (r <= split2) then
          r = r - const2
          normal_dev = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0) / &
               (((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + one)
       else
          r = r - split2
          normal_dev = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0) / &
               (((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + one)
       end if
       if (q < zero) normal_dev = - normal_dev
       return
    end if
  end subroutine standard_qnorm

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
  !
  !     Gamma distribution
  !  
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function rgamma(npar, par, rng) result(y)
    implicit none
    type(rng_t), intent(inout) :: rng
    integer, intent(in) :: npar ! must be 2
    real(dp), intent(in) :: par(npar)
    real(dp) :: y
    ! shape = nu
    ! scale = mu/nu    
    y = rgamma_default(par(2), par(1)/par(2), rng)
    return
  end function rgamma

  function rgamma_default(shape,scale,rng) result(y)
    !**********************************************************************
    !     Generates random deviates from the gamma distribution whose density is
    !     
    !          1/((scale**shape)*Gamma(shape)) * X**(shape-1) * Exp(-X/scale)
    !
    !     Arguments  (same as the rgamma in R)
    !     scale --> Location parameter of Gamma distribution  scale > 0 )
    !     shape --> Shape parameter of Gamma distribution ( shape > 0 )
    !**********************************************************************
    implicit none
    type(rng_t), intent(inout) :: rng
    real (dp), intent(in) :: shape, scale
    real (dp) :: y

    y = 1.d0 ! inicialization to avoid compiler warning
    if (scale<=0.0 .or. shape<=0.0) then
       return
    end if
    y = random_standard_gamma(shape,rng)*scale
    return    
  end function rgamma_default

  function random_standard_gamma(a,rng) result(fn_val)  
    !**********************************************************************
    !     (standard-)  gamma  distribution
    !**********************************************************************
    !               parameter  a >= 1.0
    !
    !     for details see:
    !
    !               Ahrens, J.H. and Dieter, U.
    !               Generating gamma variates by a
    !               modified rejection technique.
    !               comm. acm, 25,1 (jan. 1982), 47 - 54.
    !
    !     step numbers correspond to algorithm 'gd' in the above paper
    !                                 (straightforward implementation)
    !
    !**********************************************************************
    !               parameter  0.0 < a < 1.0
    !
    !     for details see:
    !
    !               Ahrens, J.H. and Dieter, U.
    !               Computer methods for sampling from gamma,
    !               beta, poisson and binomial distributions.
    !               computing, 12 (1974), 223 - 246.
    !
    !     (adapted implementation of algorithm 'gs' in the above paper)
    !
    !**********************************************************************
    !
    ! Modified by Taiane S. Prass, January 10, 2013, to remove the "save"
    ! variables. Those values are now storaged in the rng_t type variable
    ! so the subroutine can be used in parallel computing.
    !
    !     input: a = parameter (mean) of the standard gamma distribution
    !     rng = rng_t type variable with auxiliar parameters
    !     output: fn_val = sample from the gamma-(a)-distribution
    !
    implicit none
    type(rng_t), intent(inout) :: rng
    real (dp), intent (in) :: a
    real (dp)  :: fn_val
    !     sqrt32 is the squareroot of 32 = 5.656854249492380
    !     coefficients q(k) - for q0 = sum(q(k)*a**(-k))
    !     coefficients a(k) - for q = q0+(t*t/2)*sum(a(k)*v**k)
    !     coefficients e(k) - for exp(q)-1 = sum(e(k)*q**k)
    real(dp), parameter ::  sqrt32 = 5.656854
    real (dp), parameter :: a1 = 0.3333333
    real (dp), parameter :: a2 = -0.2500030
    real (dp), parameter :: a3 = 0.2000062
    real (dp), parameter :: a4 = -0.1662921
    real (dp), parameter :: a5 = 0.1423657
    real (dp), parameter :: a6 = -.1367177
    real (dp), parameter :: a7 = 0.1233795
    real (dp), parameter :: q1 = 0.04166669
    real (dp), parameter :: q2 =  0.02083148
    real (dp), parameter :: q3 = 0.00801191
    real (dp), parameter :: q4 =  0.00144121
    real (dp), parameter :: q5 = -0.00007388
    real (dp), parameter :: q6 =  0.00024511
    real (dp), parameter :: q7 = 0.00024240
    real (dp), parameter :: e1 = 1.
    real (dp), parameter :: e2 =  0.4999897
    real (dp), parameter :: e3 = 0.1668290
    real (dp), parameter :: e4 =  0.0407753
    real (dp), parameter :: e5 = 0.0102930
    real (dp) :: b0, e, p, q, r, t, u, v, w, x

    !   /* --- a >= 1 : GD algorithm --- */
    !    /* Step 1: Recalculations of s2, s, d if a has changed */
    if (a /= rng%aa) then
       if (a < 1.0) go to 40
       rng%aa = a
       rng%s2 = a - 0.5
       rng%s = sqrt(rng%s2)
       rng%d = sqrt32 - 12.0*rng%s
    end if

    !   /* Step 2: t = standard normal deviate,
    !               x = (s,1/2) -normal deviate. */
    !    /* immediate acceptance (i) */
    t = random_standard_norm(rng)
    x = rng%s + 0.5*t
    fn_val = x*x
    if (t >= 0.0) return

    ! /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
    u = rng_uniform(rng)
    if (rng%d*u <= t*t*t) return

    ! /* Step 4: recalculations of q0, b, si, c if necessary */
    if (a /= rng%aaa) then
       rng%aaa = a
       r = 1.0/a
       rng%q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r

       !/* Approximation depending on size of parameter a */
       !/* The constants in the expressions for b, si and c */
       !/* were established by numerical experiments */
       if (a > 3.686) then
          if (a > 13.022) then
             !/*  CASE 3:  A .GT. 13.022 */
             rng%b = 1.77
             rng%si = 0.75
             rng%c = 0.1515/rng%s
             go to 10
          end if
          !/*  CASE 2:  3.686 .LT. A .LE. 13.022 */
          rng%b = 1.654 + 0.0076*rng%s2
          rng%si = 1.68/rng%s + 0.275
          rng%c = 0.062/rng%s + 0.024
          go to 10
       end if
       !/*  CASE 1:  A .LE. 3.686 */
       rng%b = 0.463 + rng%s + 0.178*rng%s2
       rng%si = 1.235
       rng%c = 0.195/rng%s - 0.079 + 0.16*rng%s
    end if

10  continue
    !/* Step 5: no quotient test if x not positive */
    if (x > 0.0) then
       !/* Step 6: calculation of v and quotient q */
       v = t/(rng%s + rng%s)
       if (abs(v) > 0.25) then
          q = rng%q0 - rng%s*t + 0.25*t*t + (rng%s2+rng%s2)*log(1.0+v)
       else
          q = rng%q0 + 0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+ a3)*v+a2)*v+a1)*v
       end if

       !/* Step 7: quotient acceptance (q) */
       if (log(1.0-u) <= q) return
    end if

20  continue

    !/* Step 8: e = standard exponential deviate
    !*   u =  0,1 -uniform deviate
    !*   t = (b,si)-double exponential (laplace) sample */
    !/* Step    9:  rejection if t < tau(1) = -0.71874483771719 */
    e = random_standard_exponential(rng)
    u = rng_uniform(rng)
    u = u + u - 1.0
    t = rng%b + sign(rng%si*e,u)

    !/*   STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719 */
    if (t < -0.71874483771719) goto 20

    !/* Step 10:    calculation of v and quotient q */
    v = t/(rng%s + rng%s)
    if (abs(v) > 0.25) then
       q = rng%q0 - rng%s*t + 0.25*t*t + (rng%s2 + rng%s2)*log(1.0+v)
    else
       q = rng%q0 + 0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
    end if

    !/* Step 11:    hat acceptance (h) */
    !/* (if q not positive go to step 8) */
    !/* repeat .. until  `t' is accepted */
    if (q <= 0.0) go to 20

    if (q > 0.5) then
       ! if (q >= 15.0) then
       ! if (q + e - 0.5*t*t > 87.49823) go to 30
       ! if (rng%c*abs(u) > exp( q + e - 0.5*t*t)) go to 20
       ! go to 30
       ! end if
       w = exp(q) - 1.0
    else
       w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q
    end if

    !/* IF T IS REJECTED, SAMPLE AGAIN AT STEP 8 */
    if (rng%c*abs(u) > w*exp(e-0.5*t*t)) go to 20

30  continue

    x = rng%s + 0.5*t
    fn_val = x*x
    return

40  continue
    !/* GS algorithm for parameters a < 1 */
    b0 = 1.0 + 0.3678794*a

50  continue

    p = b0*rng_uniform(rng)
    if ( p < 1.0) then
       fn_val = exp(log(p)/a)
       if (random_standard_exponential(rng) < fn_val)  go to 50
       return
    end if
    fn_val = -log((b0-p)/a)
    if (random_standard_exponential(rng) < (1.0-a)*log(fn_val)) go to 50
    return
  end function random_standard_gamma

  function d_gamma(x, npar, par, give_log) result(fn)
    implicit none
    integer, intent(in) :: npar
    real(dp), intent(in) :: x, par(npar)
    logical, intent(in) :: give_log
    real(dp) :: fn      
    fn = dgamma_default(x, par(2), par(1)/par(2), give_log)
    return
  end function d_gamma

  function dgamma_default(x, shape, scale, give_log) result(fn)
    !-----------------------------------------------------
    ! Density funtion - Gamma distribution
    ! 
    ! Implemented by Taiane S. Prass
    ! April, 2020
    !-----------------------------------------------------
    !
    !  1/((scale**shape)*Gamma(shape)) * X**(shape-1) * Exp(-X/scale)
    !
    ! based on R code dgamma.c
    !
    implicit none
    real (dp), intent(in) :: x, shape, scale
    logical, intent(in) :: give_log
    real (dp) :: fn

    fn = 0.d0
    if(give_log) fn = -Huge(1.d0)

    if(x < 0) return

    if (shape == 0.0d0) then
       !      /* point mass at 0 */
       if(x == 0.0d0) fn = Huge(1.d0)
       return
    end if

    if (x == 0.0d0) then
       if(shape < 1) fn = Huge(1.d0)
       if(shape == 0.0d0) then
          if(give_log) fn = -log(scale)
          if(.not. give_log) fn = 1.d0 / scale               
       end if
       return
    end if
    !
    ! fn = -shape*log(scale) - log(gamma(shape)) + (shape - 1)*log(x) -x/scale
    ! fn = exp(fn)
    !
    if(shape < 1) then
       fn = dpois_raw(shape, x/scale, give_log)
       if(give_log) then
          !/* NB: currently *always*  shape/x > 0  if shape < 1:
          !* -- overflow to Inf happens, but underflow to 0 does NOT : */
          if(shape/x < Huge(1.d0)) then
             fn = fn + log(shape/x)
          else 
             ! /* shape/x overflows to +Inf */
             fn = fn + log(shape) - log(x)
          end if
       else
          fn = fn*shape/x
       end if
       return
    end if

    ! /* else  shape >= 1 */
    fn = dpois_raw(shape - 1.d0, x/scale, give_log)
    if(give_log) then
       fn = fn - log(scale)
    else 
       fn = fn/scale
    end if
    return 
  end function dgamma_default

  function dpois_raw(x, lambda, give_log) result(fn)
    !
    ! Adapted from R code dpois.c 
    !
    ! *    dpois_raw() computes the Poisson probability  lb^x exp(-lb) / x!.
    ! *      This does not check that x is an integer, since dgamma() may
    ! *      call this with a fractional x argument. Any necessary argument
    ! *      checks should be done in the calling function.
    ! *
    ! */
    !
    !----------------------------------------------------------------
    ! FORTRAN version by 
    ! Taiane Schaedler Prass
    ! Porto Alegre, April, 2020
    implicit none
    real(dp), intent(in) :: x, lambda
    logical, intent(in) :: give_log
    real(dp) :: fn

    fn = 0.d0
    if(give_log) fn = -Huge(1.d0)

    if(lambda == 0.0d0) then
       if(x == 0.0d0) then
          fn = 1.d0				
          if(give_log) fn = 0.d0				
       end if
       return
    end if

    ! including for the case where  x = lambda = +Inf
    if(lambda > Huge(1.d0)) return

    if(x < 0.d0)  return

    if (x <= lambda * tiny(1.d0)) then
       fn = -lambda
       if(.not. give_log) fn = exp(fn)
       return
    end if

    if (lambda < x * tiny(1.d0)) then         
       if (x > Huge(1.d0)) return !// lambda < x = +Inf
       !// else 
       fn = -lambda + x*log(lambda) - lngamma(x+1.d0)
       if(.not. give_log) fn = exp(fn)
       return
    end if

    if(give_log) then 
       fn = -0.5d0*log(M_2PI*x) -stirlerr(x)-bd0(x,lambda)
       return
    end if
    fn = exp(-stirlerr(x)-bd0(x,lambda))/sqrt(M_2PI*x)		
    return		
  end function dpois_raw

  function bd0(x, np) result(fn)
    !
    ! Adapted from R code bd0.c 
    !
    !*   Evaluates the "deviance part"
    !*   bd0(x,M) :=  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =
    !*        =  x * log(x/M) + M - x
    !*   where M = E[X] = n*p (or = lambda), for     x, M > 0
    !
    !-----------------------------------------------------------
    ! FORTRAN version by 
    ! Taiane Schaedler Prass
    ! Porto Alegre, April, 2020
    !
    implicit none
    real(dp), intent(in) :: x, np
    real(dp) :: fn
    real(dp) :: v, s, ej, s1
    integer :: j

    if(abs(x - np) < 0.1d0*(x + np)) then
       v = (x - np)/(x + np) !// might underflow to 0
       s = (x - np)*v  !/* s using v -- change by MM */

       if(abs(s) < tiny(1.d0)) then !DBL_MIN
          fn = s
          return
       end if

       ej  = 2.d0*x*v
       v  = v*v 

       do j = 1,999  
          !#/* Taylor series; 1000: no infinite loop
          !# as |v| < .1,  v^2000 is "zero" */
          ej = ej* v !// = v^(2j+1)
          s1 = s + ej/dble(2*j+1)
          if (s1 == s) then
             !/* last term was effectively 0 */
             fn = s1
             return
          end if
          s = s1
       end do
    end if

    !/* else:  | x - np |  is not too small */
    fn = x * log(x/np) + np - x
    return
  end function bd0

  function stirlerr(n) result(fn)
    !
    ! Adapted from R code stirlerr.c
    !
    !*  DESCRIPTION
    !*
    !*    Computes the log of the error term in Stirling's formula.
    !*      For n > 15, uses the series 1/12n - 1/360n^3 + ...
    !*      For n <=15, integers or half-integers, uses stored values.
    !*      For other n < 15, uses lgamma directly (don't use this to
    !*        write lgamma!)
    !*
    !/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
    !*             = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
    !*             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
    !*
    !* see also lgammacor() in ./lgammacor.c  which computes almost the same!
    !*/
    !------------------------------------------------------------------
    ! FORTRAN version by 
    ! Taiane Schaedler Prass
    ! Porto Alegre, April, 2020
    !         
    implicit none
    real(dp), intent(in) :: n
    real(dp) :: fn
    real(dp) :: nn
    real(dp), parameter :: S0 = 8.3333333333333333333d-2   !/* 1/12 */
    real(dp), parameter :: S1 = 2.77777777777777777778d-3    !/* 1/360 */
    real(dp), parameter :: S2 = 7.9365079365079365079365d-4  !/* 1/1260 */
    real(dp), parameter :: S3 = 5.95238095238095238095238d-4 !/* 1/1680 */
    real(dp), parameter :: S4 = 8.417508417508417508417508d-4!/* 1/1188 */ 
    ! /*
    ! error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
    ! */      
    real(dp), parameter :: sferr_halves(31) = (/0.0d0, & !/* n=0 - wrong, place holder only */        
         1.534264097200273452913848d-1, & !/* 0.5 */
         8.10614667953272582196702d-2, & !/* 1.0 */
         5.48141210519176538961390d-2,&  !/* 1.5 */
         4.13406959554092940938221d-2, & !/* 2.0 */
         3.316287351993628748511048d-2,& !/* 2.5 */
         2.767792568499833914878929d-2,& !/* 3.0 */
         2.374616365629749597132920d-2,& !/* 3.5 */
         2.079067210376509311152277d-2,& !/* 4.0 */
         1.848845053267318523077934d-2,& !/* 4.5 */
         1.664469118982119216319487d-2,& !/* 5.0 */
         1.513497322191737887351255d-2,& !/* 5.5 */
         1.387612882307074799874573d-2,& !/* 6.0 */
         1.281046524292022692424986d-2,& !/* 6.5 */
         1.189670994589177009505572d-2, & !/* 7.0 */
         1.110455975820691732662991d-2,& !/* 7.5 */
         1.0411265261972096497478567d-2,& !/* 8.0 */
         9.799416126158803298389475d-3,& !/* 8.5 */
         9.255462182712732917728637d-3,& !/* 9.0 */
         8.768700134139385462952823d-3,& !/* 9.5 */
         8.330563433362871256469318d-3,& !/* 10.0 */
         7.934114564314020547248100d-3,& !/* 10.5 */
         7.573675487951840794972024d-3,& !/* 11.0 */
         7.244554301320383179543912d-3,& !/* 11.5 */
         6.942840107209529865664152d-3,& !/* 12.0 */
         6.665247032707682442354394d-3,& !/* 12.5 */
         6.408994188004207068439631d-3,& !/* 13.0 */
         6.171712263039457647532867d-3,& !/* 13.5 */
         5.951370112758847735624416d-3,& !/* 14.0 */
         5.746216513010115682023589d-3,& !/* 14.5 */
         5.554733551962801371038690d-3/)  !/* 15.0 */
    real(dp), parameter :: M_LN_SQRT_2PI = 9.18938533204672741780329736406d-1   !/* log(sqrt(2*pi)) == log(2*pi)/2 */

    if (n <= 15.0d0) then
       nn = n + n
       if(ceiling(nn) == floor(nn)) then
          fn = sferr_halves(ceiling(nn) + 1)
       else !## M_LN_SQRT_2PI = ln(sqrt(2*pi)) = 0.918938..
          fn = lngamma(n + 1.d0) - (n + 0.5d0)*log(n) + n - M_LN_SQRT_2PI
       end if
       return
    end if

    nn  = n*n
    if(n > 500) then
       fn = (S0-S1/nn)/n
    else if (n > 80) then
       fn = (S0-(S1-S2/nn)/nn)/n
    else if (n > 35) then
       fn = (S0-(S1-(S2-S3/nn)/nn)/nn)/n
    else 
       !## 15 < n <= 35 :
       fn = (S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n
    end if
    return
  end function stirlerr


  function duw(y, npar, par, give_log) result(fn)
    !-----------------------------------------------------
    ! Implements the density funtion related to the Unit Weibull 
    ! distribution parameterized as in Pumi and Prass (2022)
    ! 
    ! par = c(mu, lambda, rho)
    !
    ! 0 < rho < 1,    0 < mu < 1 is a scale parameter
    ! and lambda > 0 is a shape parameter.
    ! In this parameterization, mu = rho-th quantile of the
    ! Unit-Weibull distribution with shape parameter lambda
    !
    ! Implemented by Guilherme Pumi
    ! September, 2021
    !-----------------------------------------------------
    implicit none
    integer, intent(in):: npar
    real(dp), intent(in) :: y, par(npar)
    logical, intent(in) :: give_log      
    real(dp) :: rho, mu, lambda, A,  fn

    rho = par(3)
    lambda = par(2)
    mu = par(1)
    A = log(y)/log(mu)

    fn = log(lambda) - log(y) + log(log(rho)/log(mu))    
    fn = fn + (lambda - 1)*log(A) + log(rho)*(A**lambda)
    if(.not. give_log) fn = exp(fn)
    return
  end function duw

  !-----------------------------------------------------
  ! Random generation of Unit Weibull distribution parameterized 
  ! as in Pumi and Prass (2022) using the inversion method
  !
  ! par = c(mu, lambda, rho)
  !
  ! 0 < rho <1, 0< mu < 1, lambda>0 
  !-----------------------------------------------------
  function ruw(npar, par, rng) result(fn)
    ! par = c(mu, lambda, rho)
    implicit none
    integer, intent(in) :: npar   ! must be 3
    real(dp), intent(in) :: par(npar)
    type(rng_t), intent(inout) :: rng
    real(dp) :: fn, rho, lambda, mu, u

    rho = par(3)
    lambda = par(2)
    mu = par(1)
    u = rng_uniform(rng)

    fn = log(mu)*(log(u)/log(rho))**(1/lambda)
    fn = exp(fn)

    return
  end function ruw


end module RNG_mod
