!
! Copyright (c) 2019 V.Shishkin
!

module newton_methods
  use my_prec
  use matrixopr
  use solve_methods
  implicit none
contains
  ! Calculating vector of functions values using X !
  function calc_step_vector(X) result(res)
    real(mp), dimension(:) :: X
    real(mp), allocatable, dimension(:) :: res
    integer(mp) i

    allocate(res(10)); res = 0

    res(1) = cos(X(1)*X(2)) - exp(-3.0_mp*X(3)) + X(4)*X(5)**2 - X(6) - sinh(2.0_mp*X(8))*X(9) + 2*X(10) + 2.0004339741653854440_mp
    res(2) = sin(X(1)*X(2)) + X(3)*X(9)*X(7) - exp(-X(10)+X(6))
    res(2) = res(2) + 3.0_mp*X(5)**2 - X(6)*(X(8) + 1.0_mp) + 10.886272036407019994_mp
    res(3) = x(1) - X(2) + X(3) - X(4) + X(5) - X(6) + X(7) - X(8) +X(9) - X(10) - 3.1361904761904761904_mp
    res(4) = 2*cos(-X(9) + X(4)) + X(5)/(X(3) + X(1)) - sin(X(2)**2) + cos(X(7)*X(10))**2 - X(8) - 0.1707472705022304757_mp
    res(5) = sin(X(5))+ 2.0_mp*X(8)*(X(3)+X(1)) - exp(-X(7)*(-X(10)+X(6)))
    res(5) = res(5) + 2.0_mp*cos(X(2)) - 1/(X(4) - X(9)) - 0.3685896273101277862_mp
    res(6) = exp(X(1) - X(4) - X(9)) + X(5)**2/X(8) + 0.5_mp*cos(3.0_mp*X(10)*X(2)) - X(6)*X(3) + 2.0491086016771875115_mp
    res(7) = X(2)**3*X(7) - sin(X(10)/X(5) + X(8)) + (X(1) - X(6))*cos(X(4)) + X(3) - 0.7380430076202798014_mp
    res(8) = X(5)*(X(1)-2.0_mp*X(6))**2 - 2.0_mp*sin(-X(9)+X(3)) + 1.5_mp*X(4) - exp(X(2)*X(7)+X(10)) + 3.5668321989693809040_mp
    res(9) = 7.0_mp/X(6) + exp(X(5)+X(4)) - 2.0_mp*X(2)*X(7)*X(8)*X(10) + 3.0_mp*X(9) - 3.0_mp*X(1) - 8.4394734508383257499_mp
    res(10) = X(10)*X(1) +X(9)*X(2) - X(8)*X(3) + sin(X(4)+X(5)+X(6))*X(7) - 0.78238095238095238096_mp

  end function calc_step_vector

end module newton_methods
