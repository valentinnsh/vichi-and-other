!
! Copyright (c) 2019 V.Shishkin
!

module integration
  use my_prec
  use matrixopr
  use solve_methods
  implicit none
contains

  function f(x) result(res)
    real(mp) :: res, x
    res = 0.5_mp*cos(3.0_mp*x)*exp(0.4_mp*x) + 4.0_mp*sin(3.5_mp*x)*exp(-3.0_mp*x) + 3.0_mp*x
  end function f

  ! Calculating function moments
  function calc_moments(l, r, alpha, a) result(moments)
    real(mp), dimension(1:3) :: moments
    real(mp) :: l, r, alpha, a

    moments(1) = ((r-a)**(1-alpha) - (l-a)**(1-alpha))/(1-alpha)
    moments(2) = ((r-a)**(2-alpha) - (l-a)**(2-alpha))/(2-alpha) + a*moments(1)
    moments(3) = ((r-a)**(3-alpha) - (l-a)**(3-alpha))/(3-alpha) + 2*a*moments(2) - a**2*moments(1)
  end function calc_moments


  !
  function calc_step_skf(z0,z1,alpha, h, a) result(res)
    real(mp) :: res, z0,z1,alpha,h, a
    real(mp), dimension(1:3) :: moments, Ai

    moments = calc_moments(z0, z1, alpha, a)

    Ai(1) = moments(3) - moments(2)*((z0+h)+z1)+moments(1)*z1*(z0+h)
    Ai(1) = Ai(1)/(z0+h - z0)/(z1-z0)
    Ai(2) = -(moments(3)-moments(2)*(z0+z1)+moments(1)*z0*z1)/h/(z1-z0-h)
    Ai(3) = (moments(3)-moments(2)*(z0+h+z0)+moments(1)*(z0+h)*z0)/(z1-z0-h)/(z1-z0)


    res = Ai(1)*f(z0) + Ai(2)*f(z0+h) + Ai(3)*f(z1)
  end function calc_step_skf


  ! Richards checking
  function richards(s1, s2, s3) result(res)
    real(mp) :: res, s3,s2, s1

    res = abs((s3-s2)/((s2-s1)/(s3-s2)-1))
  end function richards


  function newton_koss(a, b, alpha,l) result(res)
    real(mp), dimension(1:3) :: moments
    real(mp) :: b, alpha, a, s1, s2, s3, l, res, h
    integer(mp) :: k, i, j

  end function newton_koss


end module integration
