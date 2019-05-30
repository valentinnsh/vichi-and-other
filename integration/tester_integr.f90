!
! Copyright (c) 2019 V.Shishkin
!
program DEFAULT_NAME
  use my_prec
  use matrixopr
  use solve_methods
  use integration
  implicit none

  real(mp) :: true_res


  true_res = 8.565534243

  print *, 'Метод Ньютона-Котса дает результат :'
  print *, newton_kots(1.1_mp,2.3_mp,0.4_mp,2.0_mp)

  print *, 'Метод Гаусса дает результат :'
  print *, gauss_integration(1.1_mp,2.3_mp,0.4_mp,2.0_mp)

end program DEFAULT_NAME
