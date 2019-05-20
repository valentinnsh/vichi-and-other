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
  function calc_fun_vector(X) result(res)
    real(mp), dimension(:) :: X
    real(mp), allocatable, dimension(:) :: res

    allocate(res(10)); res = 0

    res(1) = cos(X(1)*X(2)) - exp(-3.0_mp*X(3)) + X(4)*X(5)**2 - X(6) - sinh(2.0_mp*X(8))*X(9)
    res(1) = res(1) + 2.0_mp*X(10) + 2.0004339741653854440_mp
    res(2) = sin(X(1)*X(2)) + X(3)*X(9)*X(7) - exp(-X(10)+X(6))
    res(2) = res(2) + 3.0_mp*X(5)**2 - X(6)*(X(8) + 1.0_mp) + 10.886272036407019994_mp
    res(3) = X(1) - X(2) + X(3) - X(4) + X(5) - X(6) + X(7) - X(8) +X(9) - X(10) - 3.1361904761904761904_mp
    res(4) = 2*cos(-X(9) + X(4)) + X(5)/(X(3) + X(1)) - sin(X(2)**2) + (cos(X(7)*X(10)))**2 - X(8) - 0.1707472705022304757_mp
    res(5) = sin(X(5))+ 2.0_mp*X(8)*(X(3)+X(1)) - exp(-X(7)*(-X(10)+X(6)))
    res(5) = res(5) + 2.0_mp*cos(X(2)) - 1/(X(4) - X(9)) - 0.3685896273101277862_mp
    res(6) = exp(X(1) - X(4) - X(9)) + X(5)**2/X(8) + 0.5_mp*cos(3.0_mp*X(10)*X(2)) - X(6)*X(3) + 2.0491086016771875115_mp
    res(7) = X(2)**3*X(7) - sin(X(10)/X(5) + X(8)) + (X(1) - X(6))*cos(X(4)) + X(3) - 0.7380430076202798014_mp
    res(8) = X(5)*(X(1)-2.0_mp*X(6))**2 - 2.0_mp*sin(-X(9)+X(3)) + 1.5_mp*X(4) - exp(X(2)*X(7)+X(10)) + 3.5668321989693809040_mp
    res(9) = 7.0_mp/X(6) + exp(X(5)+X(4)) - 2.0_mp*X(2)*X(7)*X(8)*X(10) + 3.0_mp*X(9) - 3.0_mp*X(1) - 8.4394734508383257499_mp
    res(10) = X(10)*X(1) +X(9)*X(2) - X(8)*X(3) + sin(X(4)+X(5)+X(6))*X(7) - 0.78238095238095238096_mp

  end function calc_fun_vector

  function calc_jacobian(x) result (j)
    real(mp), dimension(:) :: x
    real(mp), allocatable, dimension(:,:) :: j

    allocate(j(10,10)); j = 0

    !----------------------------Row 1------------!
    j(1,1) = - sin(x(1)*x(2))*x(2); j(1,2) = -sin(x(1)*x(2))*x(1); j(1,3) = 3.0_mp*exp(-3.0_mp*x(3)); j(1,4) = x(5)**2
    j(1,5) = 2.0_mp*x(4)*x(5); j(1,6) = -1.0_mp; j(1,8) = -2.0_mp*cosh(2.0_mp*x(8))*x(9); j(1,9) = -sinh(2.0_mp*x(8))
    j(1,10) = 2.0_mp
    !----------------------------Row 2------------!
    j(2,1) = cos(x(1)*x(2))*x(2); j(2,2) = cos(x(1)*x(2))*x(1); j(2,3) = x(9)*x(7); j(2,5) = 6.0_mp*x(5);
    j(2,6) = -exp(-x(10)+x(6)) - x(8) - 1.0_mp; j(2,7) = x(3)*x(9); j(2,8) = -x(6); j(2,9) = x(3)*x(7);
    j(2,10) = exp(-x(10)+x(6))
    !----------------------------Row 3------------!
    j(3,1) = 1.0_mp; j(3,2) = -1.0_mp; j(3,3) = 1.0_mp; j(3,4) = -1.0_mp; j(3,5) = 1.0_mp; j(3,6) = -1.0_mp
    j(3,7) = 1.0_mp; j(3,8) = -1.0_mp; j(3,9) = 1.0_mp; j(3,10) = -1.0_mp
    !----------------------------Row 4------------!
    j(4,1) = -x(5)/(x(3)+x(1))**2; j(4,2) = -2.0_mp*cos(x(2)**2)*x(2); j(4,3) = j(4,1); j(4,4) = -2.0_mp*sin(-x(9)+x(4))
    j(4,5) = 1.0_mp/(x(3)+x(1)); j(4,7) = -2.0_mp*cos(x(7)*x(10))*sin(x(7)*x(10))*x(10); j(4,8) = -1.0_mp;
    j(4,9) = 2.0_mp*sin(-x(9)+x(4)); j(4,10) = -2.0_mp*cos(x(7)*x(10))*sin(x(7)*x(10))*x(7)
    !----------------------------Row 5------------!
    j(5,1) = 2.0_mp*x(8); j(5,2) = -2.0_mp*sin(x(2)); j(5,3) = j(5,1); j(5,4) = 1/(-x(9)+x(4))**2
    j(5,5) = cos(x(5)); j(5,6) = x(7)*exp((-x(7))*(-x(10)+x(6))); j(5,7) = (-x(10)+x(6))*exp((-x(7))*(-x(10)+x(6)))
    j(5,8) = 2.0_mp*x(3)+2.0_mp*x(1); j(5,9) = -j(5,4); j(5,10) = -j(5,6)
    !----------------------------Row 6------------!
    j(6,1) = exp(x(1)-x(4)-x(9)); j(6,2) = -1.5_mp*sin(3.0_mp*x(2)*x(10))*x(10); j(6,3) = -x(6); j(6,4) = -j(6,1);
    j(6,5) = 2.0_mp*x(5)/x(8); j(6,6) = -x(3); j(6,8) = -(x(5)/x(8))**2; j(6,9) = j(6,4);
    j(6,10) = -1.5_mp*sin(3.0_mp*x(2)*x(10))*x(2)
    !----------------------------Row 7------------!
    j(7,1) = cos(x(4)); j(7,2) = 3.0_mp*x(2)**2*x(7); j(7,3) = 1.0_mp; j(7,4) = -(x(1)-x(6))*sin(x(4))
    j(7,5) = cos(x(10)/x(5)+x(8))*x(10)*x(5)**(-2); j(7,6) = -j(7,1); j(7,7) = x(2)**3;
    j(7,8) = -cos(x(10)/x(5)+x(8)); j(7,10) = -cos(x(10)/x(5)+x(8))/x(5)
    !----------------------------Row 8------------!
    j(8,1) = 2.0_mp*x(5)*(x(1)-2.0_mp*x(6)); j(8,2) = -x(7)*exp(x(2)*x(7)+x(10)); j(8,3) = -2.0_mp*cos(-x(9)+x(3));
    j(8,4) = 1.5_mp; j(8,5) = (x(1)-2.0_mp*x(6))**2; j(8,6) = -2.0_mp*j(8,1); j(8,7) = -x(2)*exp(x(2)*x(7)+x(10));
    j(8,9) = -j(8,3); j(8,10) = -exp(x(2)*x(7)+x(10))
    !----------------------------Row 9------------!
    j(9,1) = -3.0_mp; j(9,2) = -2.0_mp*x(7)*x(8)*x(10); j(9,4) = exp(x(5)+x(4)); j(9,5) = j(9,4);
    j(9,6) = -7.0_mp*x(6)**(-2); j(9,7) = -2.0_mp*x(2)*x(8)*x(10); j(9,8) =-2.0_mp*x(7)*x(2)*x(10);
    j(9,9) = 3.0_mp; j(9,10) = -2.0_mp*x(7)*x(8)*x(2)
    !----------------------------Row 9------------!
    j(10,1) = x(10); j(10,2) = x(9); j(10,3) = -x(8); j(10,4) = cos(x(4)+x(5)+x(6))*x(7);
    j(10,5) = j(10,4); j(10,6) = j(10,4); j(10,7) = sin(x(4)+x(5)+x(6)); j(10,8) = -x(3); j(10,9) = x(2); j(10,10) = x(1)

  end function calc_jacobian


  ! Classic Newton method !
  subroutine newton_method(X, iter_num)
    real(mp), dimension(:) :: X
    real(mp), allocatable, dimension(:,:) :: P, jac
    real(mp), allocatable, dimension(:) :: prev
    integer(mp) :: swaps
    integer :: iter_num, i, n

    iter_num = 0
    n = size(X)
    allocate(prev(n));allocate(P(n,n)); allocate(jac(n,n))
    jac = 0
    ! initial approximation !

    prev = 5

    do while(sqrt(sum((X-prev)**2)) > eps)
       prev = X

       jac = calc_jacobian(prev)

       P = 0; swaps = 0
       call decompose_LU(jac, P, swaps)
       X = prev - matmul(invert_matrix(jac, P), calc_fun_vector(prev))
       iter_num = iter_num + 1
    end do
  end subroutine newton_method

  ! Modified newton method
  subroutine modified_newton_method(X, iter_num, k)
    real(mp), dimension(:) :: X
    real(mp), allocatable, dimension(:,:) :: decm, P, jac
    real(mp), allocatable, dimension(:) :: prev
    integer(mp) :: swaps
    integer :: iter_num, i, k, j, n

    iter_num = 0
    n = size(X)
    allocate(prev(n));allocate(P(n,n)); allocate(jac(n,n))
    jac = 0
    ! initial approximation !

    prev = 5

    ! classic algorithm first k operations
    do while(iter_num < k)
       prev = X

       jac = calc_jacobian(prev)

       P = 0; swaps = 0
       call decompose_LU(jac, P, swaps)
       jac = invert_matrix(jac, P)
       X = prev - matmul(jac, calc_fun_vector(prev))
       iter_num = iter_num + 1
    end do
    ! modified method
    do while(sqrt(sum((X-prev)**2)) > eps)
       prev = X
       X = prev - matmul(jac,calc_fun_vector(prev))
       iter_num = iter_num + 1
    end do

  end subroutine modified_newton_method

  ! hibrid method
  subroutine hybrid_newton_method(X, iter_num, k)
    real(mp), dimension(:) :: X
    real(mp), allocatable, dimension(:,:) :: decm, P, jac
    real(mp), allocatable, dimension(:) :: prev
    integer(mp) :: swaps
    integer :: iter_num, i, k, j, n

    iter_num = 0
    n = size(X)
    allocate(prev(n));allocate(P(n,n)); allocate(jac(n,n))
    jac = 0
    ! initial approximation !

    prev = 5;


    do while(sqrt(sum((X-prev)**2)) > eps)
       prev = X
       jac = calc_jacobian(prev)

       P = 0; swaps = 0
       call decompose_LU(jac, P, swaps)
       jac = invert_matrix(jac, P)
       do i = 1, k
          prev = X
          X = prev - matmul(jac,calc_fun_vector(prev))
          iter_num = iter_num + 1
       end do
    end do
  end subroutine hybrid_newton_method

end module newton_methods
