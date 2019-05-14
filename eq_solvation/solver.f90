module solve_methods
  use my_prec
contains

  ! Процедура, проверяющая диагональное преобладание в матрице
  subroutine check_diagonal_dominance(A, n)
    real(mp), dimension(:,:) :: A
    real(mp) :: tmp_sum

    integer(mp) :: n, i

    do i=1,n
       if (abs(A(i,i)) < (sum(abs(A(i,1:n)))-abs(A(i,i)))) then
          write(*,*) 'В матрице А Диагонального преобладания нет'
          return
       endif
    enddo
  end subroutine check_diagonal_dominance

  !------------------Метод Якоби------------------------------!
  subroutine jakob_method(A, B, X)
    integer(mp) :: i, j, n
    real(mp) ::  norm_x ! норма разности X на к-м и (к-1)-м шаге
    real(mp), dimension(:,:) :: A
    real(mp), dimension(:) :: B, X
    real(mp), allocatable, dimension(:) :: G
    real(mp), allocatable, dimension(:,:) ::  Z

    call check_diagonal_dominance(A,n)
    n = size(A(1,:))

    allocate(Z(n,n))
    allocate(G(n))


    ! Посчитаем матрицу Z = D^(-1)*(D − A)
    Z = 0
    forall(i = 1:n, j = 1:n, i .ne. j) Z(i,j) = -A(i,j)/A(i,i)

    !Посчитаем теперь ыектор G
    forall(i = 1:n) G(i) = B(i)/A(i,i)

    ! Приближаем пока norm_x = |X(k) − X(k−1)| > eps
    X = 1
    norm_x = 1

    do while(norm_x > eps)
       norm_x = sqrt(sum((matmul(Z,X)+G-X)**2))
       X = matmul(Z,X)+G
    end do

    deallocate(G)
    deallocate(Z)
  end subroutine jakob_method

  !------------------Метод Зейделя----------------------------!
  subroutine seidel_method(A,B,X)
    integer(mp) ::  i, j, n
    real(mp) ::  norm_x ! норма разности X на к-м и (к-1)-м шаге
    real(mp), dimension(:,:) :: A
    real(mp), dimension(:) :: B, X
    real(mp), allocatable, dimension(:) :: Q, prevX
    real(mp), allocatable, dimension(:,:) ::  P

    n = size(A(1,:))

    call check_diagonal_dominance(A,n)
    allocate(Q(n))
    allocate(P(n,n))
    allocate(prevX(n))


    forall(i = 1:n, j = 1:n) P(i,j) = -A(i,j)/A(i,i)
    forall(i = 1:n) Q(i) = B(i)/A(i,i)

    norm_x = 1
    X = 0
    do while(norm_x > eps)
       prevX = X
       do i = 1,n
          X(i) = Q(i) + sum(X(1:i-1)*P(i,1:i-1)) + sum(X(i+1:n)*P(i,i+1:n))
       end do

       norm_x = sqrt(sum((X-prevX)**2))
    end do
  end subroutine seidel_method

  !------------------Метод Релаксации----------------------------!
  subroutine relax_method(A,B,n,X)
    integer(mp) ::  i, j, n
    integer(mp) ::  maxcoord
    real(mp), dimension(:,:) :: A
    real(mp), dimension(:) :: B, X
    real(mp), allocatable, dimension(:) :: Q
    real(mp), allocatable, dimension(:,:) ::  P

    call check_diagonal_dominance(A,n)

    allocate(Q(n))
    allocate(P(n,n))

    forall(i = 1:n, j = 1:n) P(i,j) = -A(i,j)/A(i,i)
    forall(i = 1:n) Q(i) = B(i)/A(i,i)

    X = 0 ! Начальое приближение !

    do while(maxval(abs(Q)) > eps)
       maxcoord = sum(maxloc(abs(Q)))
       X(maxcoord) = x(maxcoord) + Q(maxcoord)

       ! Вычисляем новые невязки !
       forall(i = 1:n) Q(i) = Q(i) + P(i,maxcoord)*Q(maxcoord)
    end do
  end subroutine relax_method

end module solve_methods
