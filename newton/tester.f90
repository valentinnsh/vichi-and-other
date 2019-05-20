!
! Copyright (c) 2019 V.Shishkin
!


program DEFAULT_NAME
  use my_prec
  use matrixopr
  use solve_methods
  use newton_methods
  implicit none

  real(mp), dimension(10) :: X
  integer :: iter_num
  integer(mp) :: id

  call newton_method(X, iter_num)

  id = 10
  open(id, file = 'result.dat')

  write(id,*)'_Classic Newton method_ '
  write(id,'(F16.8)')  X
  write(id,*) iter_num, ' iterations'

  close(id)
end program DEFAULT_NAME
