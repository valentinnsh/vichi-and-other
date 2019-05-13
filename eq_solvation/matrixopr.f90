module matrixopr
  use my_prec
contains
  ! Процедуры для ввода - вывода квадратных матриц !
  subroutine print_matrix(id, n, matrix)
    integer(mp) :: id, i, n
    real(mp), dimension(:,:) :: matrix
    do i = 1, n
       write(id, *) matrix(i, 1:n)
    end do
  end subroutine print_matrix

  subroutine read_matrix(id, n, matrix)
    integer(mp) :: id, i, n
    real(mp), dimension(:,:) :: matrix

    do i = 1, n
       read(id, *) matrix(i, 1:n)
    end do
  end subroutine read_matrix

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

end module matrixopr