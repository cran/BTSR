module functions
  integer,  parameter  :: dp = kind(1.d0)          ! double precision

contains
  function runiform(n) result(y)
    !
    ! generates uniform random variates using unif_rand() from R.
    !
    ! Implemented by Taiane S. Prass - PPGEst/UFRGS
    ! July, 2023 - Porto Alegre
    !
    implicit none
    integer, intent(in) :: n
    real (dp) :: y(n), unifrnd
    integer :: i
    call rndstart()
    do i = 1,n
       y(i) = unifrnd()
    end do
    call rndend()
    return
  end function runiform

  function mean(n, y) result(m)
    implicit none
    integer,  intent(in) :: n
    real(dp), intent(in) :: y(n)
    real(dp) :: m
    m = sum(y)/dble(n)
    return
  end function mean

  function dltestt(n, y, p) result(vout)
    implicit none
    integer,  intent(in)  :: n, p
    real(dp), intent(in)  :: y(n)
    real(dp) :: vout(2)
    real(dp) :: m, ym(n), s2
    real(dp) :: sum3(n-p), sum2, sum1
    integer  :: j, i, indicate
    real(dp) :: zi(p), zj(p), tem1(p)

    m = mean(n,y)
    ym = y - m
    s2 = sum(ym**2)/dble(n - p)
    sum3 = 0.d0
    sum2 = 0.d0
    do j = (p + 1),n
       sum1 = 0
       do i = (p + 1),n
          indicate = 0
          zi = ym((i - 1):(i - p))
          zj = ym((j - 1):(j - p))
          tem1 = 0
          where(zi <= zj) tem1 = 1
          if (product(tem1) == 1) indicate = 1
          sum1 = sum1 + ym(i) * indicate
       end do
       sum2 = sum2 + sum1**2
       sum3(j - p) = abs(sum1/sqrt(dble(n - p)))
    end do
    vout(1) = sum2/(s2 * (n - p)**2)
    vout(2) = maxval(sum3)/sqrt(s2)
    return
  end function dltestt
end module functions

subroutine dltest(n, y, B, p, vals)
  use functions
  implicit none
  integer  :: n, B, p
  real(dp) :: y(n), vals(4)
  real(dp) :: stats(2)
  real(dp) :: statmat(B,2)
  integer  :: i
  real(dp) :: ptest, m(n), u(n), ys(n)
  real(dp) :: tem(n), p1, p2

  stats = DLtestt(n, y, p)
  do i = 1,B
     ! Mammen test
     ptest = (sqrt(5.d0) + 1.d0)/(2.d0 * sqrt(5.d0))
     m = (-(sqrt(5.d0) - 1.d0)/2.d0)
     u = runiform(n)
     where(ptest < u) m = (sqrt(5.d0) + 1.d0)/2.d0
     ys = (y - mean(n,y)) * (m - mean(n,m))
     statmat(i,:) = DLtestt(n, ys, p)
  end do
  tem = 0.d0
  where(abs(statmat(:,1)) > abs(Stats(1))) tem = 1.d0
  p1 = mean(n, tem)
  tem = 0.d0
  where(abs(statmat(:,2)) > abs(Stats(2))) tem = 1.d0
  p2 = mean(n,tem)
  vals = [stats, p1, p2]
  return
end subroutine dltest
