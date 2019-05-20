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
  real(mp) :: t1,t2 !timers
  integer :: iter_num
  integer(mp) :: id

  id = 100
  open(id, file = 'result.dat', action = 'write')

  call cpu_time(t1)
  call newton_method(X, iter_num)
  call cpu_time(t2)

  write(id,*)'_Classic Newton method_ '
  write(id,'(F16.8)')  X
  write(id,*) iter_num, ' iterations ',
  write(id,'F16.8', append = 'no')t2-t1
  write(id,*) 'seconds'

  write(id,*) 'check: system(x) = '
  write(id,'(F16.8)')  calc_fun_vector(X)

  close(id)

end program DEFAULT_NAME
