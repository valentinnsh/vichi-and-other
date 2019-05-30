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

  ! Richards checking
  function richardson(s1, s2, s3) result(res)
    real(mp) :: res, s3,s2, s1

    res = abs((s3-s2)**2/(2*s2-s3-s1))
  end function richardson

  function solve_third_deg_eq(inp) result(res)
    real(mp), dimension(1:3) :: inp, res
    real(mp) :: q, r, a, b, c, phi, pi

    a = inp(3); b = inp(2); c = inp(1)

    pi = acos(-1.0)
    q = (a**2-3*b)/9
    r = (2*a**3-9*a*b+27*c)/54

    if (q**3 - r**2 <= 0) then
       print *, 'Есть комплексные корни'
       return
    end if


    phi = acos(r*q**(-1.5))/3.0

    !print *, phi
    res(1) = -2*sqrt(q)*cos(phi) - a/3
    res(2) = -2*sqrt(q)*cos(phi + 2*pi/3) - a/3
    res(3) = -2*sqrt(q)*cos(phi - 2*pi/3) - a/3

    !print *, res
  end function solve_third_deg_eq
  !---------------------------------------!
  !          NEWTON - KOTS METHOD         !
  !---------------------------------------!


  ! Calculating function moments !
  function calc_nk_moments(z0, z1, al, a) result(m)
    real(mp), dimension(0:2) :: m
    real(mp) :: z1, z0, al, a

    m(0) = ((z1-a)**(1-al) - (z0-a)**(1-al))/(1-al)
    m(1) = ((z1-a)**(2-al) - (z0-a)**(2-al))/(2-al) + a*m(0)
    m(2) = ((z1-a)**(3-al) - (z0-a)**(3-al))/(3-al) + 2*a*m(1) - a*a*m(0)
  end function calc_nk_moments


  ! Calculating quadratura coefficients !
  function calc_step_newt_kots_cqf(z0,z1,al, a) result(res)
    real(mp) :: res, z0,z1,al, a
    real(mp) :: zc ! center of [z0,z1]
    real(mp), dimension(0:2) :: m
    real(mp), dimension(1:3) :: Ai

    m = calc_nk_moments(z0, z1, al, a)
    zc = (z1+z0)/2

    Ai(1) = (m(2) - m(1)*(zc+z1) + m(0)*zc*z1)/(zc-z0)/(z1-z0)
    Ai(2) = -(m(2)-m(1)*(z1+z0)+m(0)*z1*z0)/(zc-z0)/(z1-zc)
    Ai(3) = (m(2) - m(1)*(zc+z0) + m(0)*zc*z0)/(z1-zc)/(z1-z0)
    res = Ai(1)*f(z0) + Ai(2)*f(zc) + Ai(3)*f(z1)
  end function calc_step_newt_kots_cqf

  function newton_kots(a, b, al, l) result(res)
    real(mp) :: b, al, a, s1, s2, s3, l, res, h
    integer(mp) :: k, i, j

    h = (b-a)
    s1 = calc_step_newt_kots_cqf(a,b,al,a)

    h = h/2
    s2 = calc_step_newt_kots_cqf(a, a+h, al, a) + calc_step_newt_kots_cqf(a+h,b,al,a)

    h = h/2; i = 0; s3 = 0
    do while((a+h*i) < b)
       i = i + 1
       s3 = s3 + calc_step_newt_kots_cqf(a+(i-1)*h, a+i*h,al,a)
    end do

    do while(richardson(s1,s2,s3) > eps)
       s1 = s2; s2 = s3;

       h = h/2; i = 0; s3 = 0
       do while((a+h*i) < b)
          i = i + 1
          s3 = s3 + calc_step_newt_kots_cqf(a+(i-1)*h, a+i*h,al,a)
       end do
    end do

    res = s3
    !print *, richards(s1,s2,s3)
    !print *, s1, s2, s3
    !res = calc_quadr_coef(a,b,al, a)
  end function newton_kots


  !---------------------------------------!
  !           GAUSS   METHOD              !
  !---------------------------------------!

  function calc_gauss_moments(z0, z1, al, a) result(m)
    real(mp), dimension(0:5) :: m
    real(mp) :: z1, z0, al, a

    m(0) = ((z1-a)**(1-al)-(z0-a)**(1-al))/(1-al)
    m(1) = ((z1-a)**(2-al)-(z0-a)**(2-al))/(2-al)+a*m(0)
    m(2) = ((z1-a)**(3-al)-(z0-a)**(3-al))/(3-al)+2*a*m(1)-a*a*m(0)
    m(3) = ((z1-a)**(4-al)-(z0-a)**(4-al))/(4-al)+3*a*m(2)-3*a*a*m(1)+(a**3)*m(0)
    m(4) = ((z1-a)**(5-al)-(z0-a)**(5-al))/(5-al)+4*a*m(3)-6*a*a*m(2)+4*(a**3)*m(1)-(a**4)*m(0)
    m(5) = ((z1-a)**(6-al)-(z0-a)**(6-al))/(6-al)+5*a*m(4)-10*a*a*m(3)+10*(a**3)*m(2)- 5*(a**4)*m(1)+(a**5)*m(0)
  end function calc_gauss_moments

  function calc_gauss_moments_old(z0, z1, al, a) result(m)
    real(mp), dimension(0:5) :: m
    real(mp) :: z1, z0, al, a

    m(0) = (z1-a)**(1-al)/(1-al) - (z0-a)**(1-al)/(1-al)

    m(1) = ((z1-a)**(1-al)*(a-z1*(al-1)) - (z0-a)**(1-al)*(a-z0*(al-1)))/(1-al)/(2-al)

    m(2) =         (z0-a)**(1-al)*(2*a**2+2*a*z0*(1-al)+z0**2*(1-al)*(2-al))
    m(2) = (m(2) - (z1-a)**(1-al)*(2*a**2+2*a*z1*(1-al)+z1**2*(1-al)*(2-al)))/(al-3)/(al-2)/(al-1)

    m(3) =        (z1-a)**(1-al)*(6*a**3-6*a**2*z1*(al-1)+3*a*z1**2*(1-al)*(2-al)+z1**3*(1-al)*(2-al)*(3-al))
    m(3) = m(3) - (z0-a)**(1-al)*(6*a**3-6*a**2*z0*(al-1)+3*a*z0**2*(1-al)*(2-al)+z0**3*(1-al)*(2-al)*(3-al))
    m(3) = m(3)/(al-4)/(al-3)/(al-2)/(al-1)

    m(4) = -(((z1-a)**(1-al)*(24*a**4-24*a**3*z1*(al-1)+12*a**2*z1**2*(al**2-3*al+2) &
         -4*a*z1**3*(al**3-6*al**2+11*al-6)+z1**4*(al**4-10*al**3+35*al**2-50*al+24)))) &
         + (((z0-a)**(1-al)*(24*a**4-24*a**3*z0*(al-1)+12*a**2*z0**2*(al**2-3*al+2) &
         -4*a*z0**3*(al**3-6*al**2+11*al-6)+z0**4*(al**4-10*al**3+35*al**2-50*al+24))))
    m(4) = m(4)/(al-5)/(al-4)/(al-3)/(al-2)/(al-1)

    m(5) = (((-a+z1)**(1-al)*(120*a**5-120*a**4*(-1+al)*z1+60*a**3*(2-3*al+al**2)*z1**2-20*a**2*(-6+11*al-6*al**2+al**3)*z1**3 &
         +5*a*(24-50*al+35*al**2-10*al**3+al**4)*z1**4-(-120 + 274*al - 225*al**2 + 85*al**3 - 15*al**4 + al**5)*z1**5))  &
         - ((-a+z0)**(1-al)*(120*a**5-120*a**4*(-1+al)*z0+60*a**3*(2-3*al+al**2)*z0**2-20*a**2*(-6+11*al-6*al**2+al**3)*z0**3 &
         +5*a*(24-50*al+35*al**2-10*al**3+al**4)*z0**4-(-120 + 274*al - 225*al**2 + 85*al**3 - 15*al**4 + al**5)*z0**5))) &
         /((-6 + al)/(-5 + al)/(-4 + al)/(-3 + al)/(-2 + al)/(-1 + al))
  end function calc_gauss_moments_old

  function calc_step_gauss_cqf(z0,z1,al, a) result(res)
    real(mp) :: res, z0,z1, al, a
    real(mp) :: r1, r0 ! коэффициенты узлового многочлена
    real(mp) :: t1, t0
    real(mp) :: zc ! center of [z0,z1]
    real(mp), dimension(0:5) :: m
    real(mp), dimension(1:3) :: Ai, t
    real(mp), dimension(1:3) :: B, c
    real(mp), dimension(3,3) :: S, P
    integer(mp) :: swaps
    integer :: i
    S = 0; P = 0;

    zc = (z1-z0)/2
    m = calc_gauss_moments(z0,z1,al,a)
    !print *, '--------------------------------------------'
    !print *, m
    !print *, '--------------------------------------------'

    S(1,1) = m(0); S(1,2) = m(1); S(1,3) = m(2); B(1) = -m(3)
    S(2,1) = m(1); S(2,2) = m(2); S(2,3) = m(3); B(2) = -m(4)
    S(3,1) = m(2); S(3,2) = m(3); S(3,3) = m(4); B(3) = -m(5)

    call decompose_LU(S, P, swaps)
    call solve_eq_sys(S, P, B, c)
    !print *, c
    ! посчитаем, используя тригонометрическую формулу виетта

    t = solve_third_deg_eq(c)

    !print *, t
    ! solve final system
    S(1,1) = 1;       S(1,2) = 1;       S(1,3) = 1;       B(1) = m(0)
    S(2,1) = t(1);    S(2,2) = t(2);    S(2,3) = t(3);    B(2) = m(1)
    S(3,1) = t(1)**2; S(3,2) = t(2)**2; S(3,3) = t(3)**2; B(3) = m(2)

    call decompose_LU(S, P, swaps)
    call solve_eq_sys(S, P, B, Ai)

    res = Ai(1)*f(t(1)) + Ai(2)*f(t(2)) + Ai(3)*f(t(3))
  end function calc_step_gauss_cqf

  function gauss_integration(a, b, al, l) result(res)
    real(mp) :: b, al, a, s1, s2, s3, l, res, h
    integer(mp) :: k, i, j

    h = (b-a)
    s1 = calc_step_gauss_cqf(a,b,al,a)

    !print *, s1
    !print *, calc_step_gauss_cqf(a+h/2,b,al,a)

    h = h/2
    s2 = calc_step_gauss_cqf(a, a+h, al, a) + calc_step_gauss_cqf(a+h,b,al,a)

    !print *, s1
    h = h/2; i = 0; s3 = 0
    do while((a+h*i) < b)
       i = i + 1
       s3 = s3 + calc_step_gauss_cqf(a+(i-1)*h, a+i*h,al,a)
    end do

    !print *, s1, s2, s3
    do while(richardson(s1,s2,s3) > eps)
       s1 = s2; s2 = s3;

       h = h/2; i = 0; s3 = 0
       do while((a+h*i) < b)
          i = i + 1
          s3 = s3 + calc_step_gauss_cqf(a+(i-1)*h, a+i*h,al,a)
          print *, calc_step_gauss_cqf(a+(i-1)*h, a+i*h,al,a)
       end do

       print *, s3
    end do

    res = s3
    !print *, richards(s1,s2,s3)
    !print *, s1, s2, s3
    !res = calc_quadr_coef(a,b,al, a)
  end function gauss_integration
end module integration
