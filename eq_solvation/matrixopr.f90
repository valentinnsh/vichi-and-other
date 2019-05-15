module matrixopr
  use my_prec
contains
  ! Процедуры для ввода - вывода квадратных матриц !
  subroutine print_matrix(id, n, matrix)
    integer(mp) :: id, i, n
    real(mp), dimension(:,:) :: matrix
    do i = 1, n
       do j = 1,n
          write(id, '(F16.8)', advance="no") matrix(i, j)
       end do
       write(id,*) ''
    end do
  end subroutine print_matrix

  subroutine read_matrix(id, n, matrix)
    integer(mp) :: id, i, n
    real(mp), dimension(:,:) :: matrix

    do i = 1, n
       read(id, *) matrix(i, 1:n)
    end do
  end subroutine read_matrix

  function matr_norm(M) result(res)
    integer(mp) :: id, i, n
    real(mp), dimension(:,:) :: M
    real(mp) :: tmp, res
    n = size(M(1,:))

    do j = 1,n
       tmp = 0
       do i = 1,n
          tmp = tmp + abs(M(i,j))
       end do
       if(tmp > res) res = tmp
    end do
  end function matr_norm

  ! Процедура производит LU  разложение матрицы А !
  ! P - матрица перестановок, для решения СЛАУ    !
  ! swaps - значение четности числа перестановок  !
  subroutine decompose_LU(decA, P, swaps)
    real(mp), dimension(:,:) :: decA, P
    integer(mp) :: maxel_loc(1)
    real(mp) :: tmp
    integer(mp) :: swaps, i, j, n

    swaps = 1
    n = size(decA(1,:))
    P = 0
    forall(i = 1:n) P(i,i) = 1

    do j = 1,n

       ! Выбор ведущего элемента в j-м столбце !

       if(abs(maxval(decA(j:n,j))) < abs(minval(decA(j:n,j)))) then
          maxel_loc = minloc(decA(j:n,j)) + j-1
       else
          maxel_loc = maxloc(decA(j:n,j)) + j-1
       end if

       if(maxel_loc(1) .ne. j) then
          ! меняем местами j-ю и maxel_loc строчки !
          swaps = swaps*(-1)
          do i = 1,n
             tmp =  decA(j,i)
             decA(j,i) = decA(maxel_loc(1), i)
             decA(maxel_loc(1), i) = tmp

             tmp =  P(j,i)
             P(j,i) = P(maxel_loc(1), i)
             P(maxel_loc(1), i) = tmp
          end do
       end if

       do i = j+1, n
          decA(i,j) = decA(i,j)/decA(j,j)
          do k = j+1,n
             decA(i,k) = decA(i,k) - decA(i,j)*decA(j,k)
          end do
       end do
    end do
  end subroutine decompose_LU

  ! Функция, считающая определитель с помощью LUP разложения !
  function determinant(A) result(det)
    real(mp), dimension(:,:) :: A
    real(mp), allocatable, dimension(:,:) :: decA, P
    real(mp) :: det
    integer(mp) :: swaps, i, j, n

    n = size(A(1,:))

    allocate(decA(n,n))
    allocate(P(n,n))

    decA = A; P = 0;
    call decompose_LU(decA, P, swaps)

    det = swaps
    do i = 1,n
       det = det*decA(i,i)
    end do

  end function determinant

  subroutine solve_eq_sys(decA, P, B, X)
    real(mp), dimension(:,:) :: decA, P
    real(mp), dimension(:) :: B, X
    real(mp), allocatable, dimension(:) :: decB
    real(mp) :: det
    integer(mp) :: swaps, i, j, n

    n = size(decA(1,:))

    allocate(decB(n))
    decB = matmul(P,B)

    ! первый шаг - прямая подстановка !
    do i = 1, n
       X(i) = decB(i)
       do j = 1, i-1
          X(i) = X(i) - X(j)*decA(i,j)
       end do
    end do

    ! обратная подстановка !
    do i = n, 1, -1
       do j = n, i+1, -1
          X(i) = X(i) - X(j)*decA(i,j)
       end do
       X(i) = X(i)/decA(i,i)
    end do
  end subroutine solve_eq_sys

  ! Функция обращения матрицы. На вход подаются результаты LU разложения !
  function invert_matrix(decA, P) result(inv)
    real(mp), dimension(:,:) :: decA, P
    real(mp), allocatable,dimension(:,:) :: inv
    real(mp), allocatable, dimension(:) :: B, X
    real(mp), allocatable, dimension(:) :: decB
    integer(mp) ::  i, j, n

    n = size(decA(1,:))
    allocate(B(n)); allocate(X(n));
    allocate(inv(n,n))
    do i = 1,n
       B = 0;
       B(i) = 1

       call solve_eq_sys(decA,P,B,X)
       do j = 1,n
          inv(j,i) = X(j)
       end do
    end do
  end function invert_matrix

  ! Нахлждение числа обусловленности матрицы !
  function condition_number(A,invA) result(cond)
    real(mp), dimension(:,:) :: A, invA
    real(mp) :: cond, tmpQ, tmpI, normA, normInv
    integer(mp) ::  i, j, n

    n = size(A(1,:))
    cond = 0; normInv = 0; normA = 0;

    do i = 1,n
       tmpA = 0
       tmpI = 0
       do j = 1,n
          tmpA = tmpA + abs(A(i,j))**2
          tmpI = tmpI + abs(invA(i,j))**2
       end do
       normA = normA + sqrt(tmpA)
       normInv = normInv + sqrt(tmpI)
    end do

    cond = normInv*normA
  end function condition_number


end module matrixopr
