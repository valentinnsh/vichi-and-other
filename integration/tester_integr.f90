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


  true_res = 8.56553

  print *, newton_koss(1.1_mp,2.3_mp,0.4_mp,2.0_mp)

end program DEFAULT_NAME
