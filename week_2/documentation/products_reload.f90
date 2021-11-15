program matrix_products
    implicit none

    interface
        subroutine mat_mul_1(mat_a,mat_b,mat_c,debug)
            real*8, intent(out):: mat_c(:,:)
            real*8, intent(in):: mat_a(:,:), mat_b(:,:)
            logical:: debug
        end subroutine
    end interface
    
    interface
        subroutine mat_mul_2(mat_a,mat_b,mat_c,debug)
            real*8, intent(out):: mat_c(:,:)
            real*8, intent(in):: mat_a(:,:), mat_b(:,:)
            logical:: debug
        end subroutine
    end interface
    
    real*8, allocatable :: A(:,:), B(:,:), C(:,:), D(:,:), E(:,:)
    integer*4 :: iterations, min_N, max_N, ii, jj, kk = 1, errors, slice_size = 5
    real*4 :: start, finish
    real*4, allocatable :: timing(:,:)
    logical :: write_every_n = .False., debug = .False.
    
    write(*,*) "Enter the timing iterations and the range of rows the test matrices should have"
    read(*,*) iterations, min_N, max_N

    ! Validating that at least one timing iteration is performed.
    ! By construction, this also validates that the number of iterations is a positive integer.
    if (iterations .lt. 1) stop "The number of iterations should be at least one"
    
    ! Validating that the minimum number of rows is at least 2.
    if (min_N .lt. 2) stop "The minimum number of rows for the square matrices must be at least 2"
    ! Validating that the maximum number of rows is higher than the minimum number of rows.
    ! Applied sequentually with the validation above, this guarantees max_N is a positive integer.
    if (max_N .le. min_N) stop "The maximum number of rows for the square matrices must be higher than the minimum number of rows."
    
    ! Allocate the timing array.
    allocate(timing(min_N:max_N,0:3),stat=errors)
    ! Validate the allocation was successful.
    ! Most likely cause for unsuccessful allocation is lack of memory.
    ! A very wide gap between min_N and max_N may lead to this behaviour
    if (errors.ne.0) then
        write(*,*) "Unable to allocate array for timing."
        write(*,*) "Writing to file performances every_slice size steps."
        write_every_n = .True.
    else if (write_every_n .eqv. .False.) then
        kk = min_N
    end if

    ! Performing outside the earlier control in case write_every_n is initialized as .True. for debugging
    if (write_every_n) then
        deallocate(timing)
        allocate(timing(1:slice_size,0:3))
    end if

    ! Finally, validate the maximum matrix size is not too big for allocation
    ! Since we're allocating five matrices for the timing, we take as proxy a (5*max_N)x(max_N) array
    allocate(A(1:5*max_N,1:max_N), stat = errors)
    if (errors.ne.0) stop "ERROR: Allocation issues. Try a smaller upper bound for the number of rows."
    deallocate(A)

    ! Open the file to store the performance data.
    open(15,file = "performance.dat")

    ! Generating two square random matrices, from min_N x min_N to max_N x max_N
    do ii = min_N, max_N
        ! "input size"
        timing(kk,0) = ii
    
        ! Allocate all the matrices using the new size
        allocate(A(1:ii,1:ii))
        allocate(B(1:ii,1:ii))
        allocate(D(1:ii,1:ii))
        allocate(E(1:ii,1:ii))
        allocate(C(1:ii,1:ii))
    
        ! Populate the arrays with uniformly distributed values
        call random_number(A)
        call random_number(B)
    
    
        !----------------------
        !    Timing products
        !----------------------
        ! Measuring performance of the first product definition
        call cpu_time(start)
        do jj = 1, iterations
            call mat_mul_1(A,B,C,debug)
        end do
        call cpu_time(finish)
        timing(kk,1) = (finish-start)/iterations
    
        ! Measuring performance of the second product definition
        ! Transposing and storing before entering the timing loop
        D = transpose(B)
        E = transpose(A)
        call cpu_time(start)
        do jj = 1, iterations
            call mat_mul_2(D,E,C,debug) ! This gives us C transpose
        end do
        call cpu_time(finish)
        timing(kk,2) = (finish-start)/iterations
    
        ! Measuring performance of built-in product definition
        call cpu_time(start)
        do jj = 1, iterations
            C = matmul(A,B)
        end do
        call cpu_time(finish)
        timing(kk,3) = (finish-start)/iterations
    
        ! Release the variables for the incoming allocation
        deallocate(A)
        deallocate(B)
        deallocate(D)
        deallocate(E)
        deallocate(C)

        ! Write to file when slice_size timing rows have been reached
        ! or when we have reached the maximum number of matrix's rows.
        if (write_every_n .and. (kk .eq. slice_size .or. ii .eq. max_N)) then
            ! Write up until the collected value.
            do jj = 1, kk
                write(15,*) timing(jj,0:3)
            end do
            ! set kk to zero, so that right after the control it is set to one.
            kk = 0
        end if
        ! Increase kk, usually the same as taking ii, but handled for sliced timings.
        kk = kk + 1
    end do
    
    ! Write the performance measurements
    if (write_every_n .eqv. .False.) then
        do ii = min_N, max_N
            write(15,*) timing(ii,0:3)
        end do
    end if
    close(15)
    write(*,*) "Averaged timings of operations successfully computed and written to file."
    
end program matrix_products
    
!--------------------------------
! Subroutines for matrix products
!--------------------------------

subroutine mat_mul_1(mat_a,mat_b,mat_c, debug)
    implicit none
    real*8, intent(in):: mat_a(:,:), mat_b(:,:)
    real*8, intent(out) :: mat_c(:,:)
    integer*4:: M,N,O,P,ii,jj,kk
    real*8:: cumulative
    integer*4:: shape_a(2), shape_b(2)
    logical:: debug
    
    shape_a = shape(mat_a)
    shape_b = shape(mat_b)
    M = shape_a(1)
    N = shape_a(2)
    O = shape_b(1)
    P = shape_b(2)

    ! Validate the matrices shapes are compatible
    if (debug .and. (N .ne. O)) stop "Input matrices sizes are not compatible for multiplication"

    do ii = 1,M
        do jj = 1, P
            cumulative = 0.d0
            do kk = 1, N
                cumulative = cumulative + mat_a(ii,kk)*mat_b(kk,jj)
            end do
            mat_c(ii,jj) = cumulative
        end do
    end do
    
end subroutine
        
    
subroutine mat_mul_2(mat_a,mat_b, mat_c, debug)
    implicit none
    real*8, intent(in):: mat_a(:,:), mat_b(:,:)
    real*8, intent(out) :: mat_c(:,:)
    integer*4:: M,N,O,P,ii,jj,kk
    real*8:: cumulative
    integer*4:: shape_a(2), shape_b(2)
    logical:: debug

    shape_a = shape(mat_a)
    shape_b = shape(mat_b)
    M = shape_a(1)
    N = shape_a(2)
    O = shape_b(1)
    P = shape_b(2)

    ! Validate the matrices shapes are compatible
    if (debug .and. (N .ne. O)) stop "Input matrices sizes are not compatible for multiplication"

    do jj = 1,P
        do ii = 1, M
            cumulative = 0.d0
            do kk = 1, N
                cumulative = cumulative + mat_a(ii,kk)*mat_b(kk,jj)
            end do
            mat_c(ii,jj) = cumulative
        end do
    end do
    
end subroutine
        
    
    