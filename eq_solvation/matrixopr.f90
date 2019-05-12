module matrixopr
  use my_pec
contains
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
end module matrixopr
