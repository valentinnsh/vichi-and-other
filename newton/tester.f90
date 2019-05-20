!
! Copyright (c) 2019 V.Shishkin
!


program DEFAULT_NAME
  use my_prec
  use matrixopr
  use solve_methods
  use newton_methods
  implicit none

  real(mp), dimension(10) :: X, X1
  real(mp) :: t1,t2 !timers
  integer :: iter_num, i, k
  integer(mp) :: id

  X  = (/0.5_mp, 0.5_mp, 1.5_mp, -1.0_mp, -0.5_mp, 1.5_mp, 0.5_mp, -0.5_mp, 1.5_mp, -1.5_mp/)
  id = 100
  open(id, file = 'result.dat', action = 'write')

  call cpu_time(t1)
  call newton_method(X, iter_num)
  call cpu_time(t2)

  write(id,*)'Classic Newton method: '
  write(id,'(F16.8)')  X
  write(id,*) iter_num, ' iterations '
  write(id,'(F16.8)', advance = "no") t2-t1
  write(id,*) 'seconds'

  write(id,*) 'check: system(x) = '
  write(id,'(F16.8)')  calc_fun_vector(X)

  X  = (/0.5_mp, 0.5_mp, 1.5_mp, -1.0_mp, -0.5_mp, 1.5_mp, 0.5_mp, -0.5_mp, 1.5_mp, -1.5_mp/)
  write(id,*) 'Modified method:'

  ! Optimal k for modified is 3 !
  !do i = 1,5
  k = 3
  call cpu_time(t1)
  call modified_newton_method(X, iter_num,k)
  call cpu_time(t2)

  write(id,*) 'check: system(x) = '
  write(id,'(F16.8)')  calc_fun_vector(X)

  write(id,*) iter_num, 'iterations'
  write(id,*) 'k = ', k
  write(id, '(F16.8)', advance = "no") t2 - t1
  write(id,*) 'seconds'
  !end do

  X  = (/0.5_mp, 0.5_mp, 1.5_mp, -1.0_mp, -0.5_mp, 1.5_mp, 0.5_mp, -0.5_mp, 1.5_mp, -1.5_mp/)
  ! Optimal k for hybrid is 2 !
  k = 2
  write(id,*) 'Hybrid method: '
  !do i = 1,5

  call cpu_time(t1)
  call hybrid_newton_method(X, iter_num,k)
  call cpu_time(t2)

  write(id,*) 'check: system(x) = '
  write(id,'(F16.8)')  calc_fun_vector(X)

  write(id,*) iter_num, 'iterations'
  write(id,*) 'k = ', k
  write(id, '(F16.8)', advance = 'no') t2 - t1
  write(id,*) 'seconds'
  !end do

  X1  = (/0.5_mp, 0.5_mp, 1.5_mp, -1.0_mp, -0.2_mp, 1.5_mp, 0.5_mp, -0.5_mp, 1.5_mp, -1.5_mp/)

  write(id,*) '---------------SECOND STARTING APPROXIMATION----------------------'
  call cpu_time(t1)
  call newton_method(X1, iter_num)
  call cpu_time(t2)

  write(id,*)'Classic Newton method: '
  write(id,'(F16.8)')  X1
  write(id,*) 'check: system(x) = '
  write(id,'(F16.8)')  calc_fun_vector(X)
  write(id,*) iter_num, ' iterations '
  write(id,'(F16.8)', advance = "no") t2-t1
  write(id,*) 'seconds'


  X1  = (/0.5_mp, 0.5_mp, 1.5_mp, -1.0_mp, -0.2_mp, 1.5_mp, 0.5_mp, -0.5_mp, 1.5_mp, -1.5_mp/)
  write(id,*) 'Modified method:'

  ! Optimal k for modified is > 8 !
  k = 8
  !do i = 1,5
  call cpu_time(t1)
  call modified_newton_method(X1, iter_num,k)
  call cpu_time(t2)

  write(id,'(F16.8)')  X1

  write(id,*) iter_num, 'iterations'
  write(id,*) 'k = ', k
  write(id, '(F16.8)', advance = "no") t2 - t1
  write(id,*) 'seconds'
  !end do


  X1  = (/0.5_mp, 0.5_mp, 1.5_mp, -1.0_mp, -0.2_mp, 1.5_mp, 0.5_mp, -0.5_mp, 1.5_mp, -1.5_mp/)
  ! Optimal k for hybrid is > 8!
  k = 9
  write(id,*) 'Hybrid method: '
  !do i = 1,5
  call cpu_time(t1)
  call hybrid_newton_method(X1, iter_num,k)
  call cpu_time(t2)

  write(id,'(F16.8)')  X1

  write(id,*) iter_num, 'iterations'
  write(id,*) 'k = ', k
  write(id, '(F16.8)', advance = 'no') t2 - t1
  write(id,*) 'seconds'
  !end do

  close(id)

end program DEFAULT_NAME
