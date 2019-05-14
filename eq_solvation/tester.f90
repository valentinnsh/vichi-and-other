!
! Copyright (c) 2019 V.Shishkin
!


program Task_1
  use my_prec
  use matrixopr
  use solve_methods
  use makedata
  implicit none

  real(mp), allocatable, dimension(:,:) :: A, decA, P, L, U
  real(mp), allocatable, dimension(:) :: B,X
  real(mp) :: tmp
  integer(mp) ::  i, j, n, id, swaps, size, q

  id = 100
  size = 5
  call make_data(size)
  open(id, file='data.dat')
  read(id,'(2x,I6)') n

  allocate(B(n)); allocate(X(n))
  allocate(A(n,n)); allocate(decA(n,n))
  allocate(P(n,n)); allocate(L(n,n)); allocate(U(n,n))
  L = 0; U = 0; decA = 0; P = 0; A = 0; X = 0;

  call read_matrix(id, n, A)
  do i=1,n
     read(id,*) B(i)
  end do
  close(id)

  decA = A
  call decompose_LU(decA, P, swaps)

  forall(i = 1:n) L(i,i) = 1
  do i = 1,n
     do j = 1,i-1
        L(i,j) = decA(i,j)
     end do
     do j = i,n
        U(i,j) = decA(i,j)
     end do
  end do




  id = 200
  open(id, file = 'result.dat')
  write(id,*) 'check LU - PA = :'
  call print_matrix(id,n,matmul(L,U)-matmul(P,A))
  write(id,*) 'det(A) = ', determinant(A)

  call solve_eq_sys(decA, P, B,X)

  write(id,*) 'Ax-b = ', matmul(A,X) - B
  !call print_matrix(id,n,P)

  write(id,*) 'A**-1*A = '
  call print_matrix(id,n, matmul(invert_matrix(decA,P),A))

  write(id, *) 'A*A**-1 = '
  call print_matrix(id,n, matmul(A,invert_matrix(decA,P)))

  write(id, *) 'condition number = ||A||*||A**-1|| = ', condition_number(A,invert_matrix(decA,P))

  close(id)

  size = 5; q = 4
  call make_diag_data(size,q)

  open(id, file='data_diag.dat')
  read(id,'(2x,I6)') n

  call read_matrix(id, n, A)
  do i=1,n
     read(id,*) B(i)
  end do
  close(id)
  X = 0



end program Task_1
