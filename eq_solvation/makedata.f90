!
! Copyright (c) 2019 V.Shishkin
!

module makedata
  use matrixopr
  use my_prec
contains
  subroutine make_data(n)
    integer(mp) i, j, id, n
    real(mp) tmp
    real(mp), allocatable, dimension(:,:) :: matrix
    real(mp), allocatable, dimension(:) :: vector
    id = 10
    open(id, file = 'data.dat')
    write(id,"('# ',I6)") n
    allocate(matrix(n,n))
    allocate(vector(n))
    call RANDOM_NUMBER(matrix)
    matrix = matrix*100
    call print_matrix(id, n, matrix)

    call RANDOM_NUMBER(vector)
    vector = vector*100
    do i = 1, n
       write(id, '(F16.8)') vector(i)
    end do
    close(id)

    deallocate(matrix)
  end subroutine make_data
end module makedata
