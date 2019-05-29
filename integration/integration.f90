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

  ! Calculating function moments !
  function calc_moments(z0, z1, al, a) result(m)
    real(mp), dimension(0:2) :: m
    real(mp) :: z1, z0, al, a

    m(0) = ((z1-a)**(1-al) - (z0-a)**(1-al))/(1-al)
    m(1) = ((z1-a)**(2-al) - (z0-a)**(2-al))/(2-al) + a*m(0)
    m(2) = ((z1-a)**(3-al) - (z0-a)**(3-al))/(3-al) + 2*a*m(1) - a*a*m(0)
  end function calc_moments


  ! Calculating quadratura coefficients !
  function calc_step_cqf(z0,z1,al, a) result(res)
    real(mp) :: res, z0,z1,al, a
    real(mp) :: zc ! center of [z0,z1]
    real(mp), dimension(0:2) :: m
    real(mp), dimension(1:3) :: Ai

    m = calc_moments(z0, z1, al, a)
    zc = (z1+z0)/2

    Ai(1) = (m(2) - m(1)*(zc+z1) + m(0)*zc*z1)/(zc-z0)/(z1-z0)
    Ai(2) = -(m(2)-m(1)*(z1+z0)+m(0)*z1*z0)/(zc-z0)/(z1-zc)
    Ai(3) = (m(2) - m(1)*(zc+z0) + m(0)*zc*z0)/(z1-zc)/(z1-z0)
    res = Ai(1)*f(z0) + Ai(2)*f(zc) + Ai(3)*f(z1)
  end function calc_step_cqf


  ! Richards checking
  function richardson(s1, s2, s3) result(res)
    real(mp) :: res, s3,s2, s1

    res = abs((s3-s2)**2/(2*s2-s3-s1))
  end function richardson


  function newton_kots(a, b, al, l) result(res)
    real(mp), dimension(1:3) :: moments
    real(mp) :: b, al, a, s1, s2, s3, l, res, h
    integer(mp) :: k, i, j

    h = (b-a)/2
    s1 = calc_step_cqf(a,b,al,a)

    h = h/2
    s2 = calc_step_cqf(a, a+h, al, a) + calc_step_cqf(a+h,b,al,a)

    h = h/2; i = 0; s3 = 0
    do while((a+h*i) < b)
       i = i + 1
       s3 = s3 + calc_step_cqf(a+(i-1)*h, a+i*h,al,a)
    end do

    do while(richardson(s1,s2,s3) > eps)
       s1 = s2; s2 = s3;

       h = h/2; i = 0; s3 = 0
       do while((a+h*i) < b)
          i = i + 1
          s3 = s3 + calc_step_cqf(a+(i-1)*h, a+i*h,al,a)
       end do
    end do

    res = s3
    !print *, richards(s1,s2,s3)
    !print *, s1, s2, s3
    !res = calc_quadr_coef(a,b,al, a)
  end function newton_kots


end module integration
