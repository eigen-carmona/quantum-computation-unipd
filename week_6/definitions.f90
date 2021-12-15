module definitions
    implicit none
        type complex_matrix
            integer*4 :: dimensions(1:2)
            complex*16, allocatable :: entries(:,:)
            complex*16 :: determinant, trace
        end type
    
    contains
        function enter_mat(entries) result(new_mat)
            implicit none
            complex*16, intent(in) :: entries(:,:)
            type(complex_matrix):: new_mat
    
            new_mat%dimensions = shape(entries)
            new_mat%entries = entries
            new_mat%trace = mat_trace(entries)
            ! This function is intended for the initialization of the matrix
        end function enter_mat
    
        function mat_trace(entries) result(trace)
            implicit none
            complex*16, intent(in) :: entries(:,:)
            complex*16 :: trace
            integer*4 :: ii, NN
            integer*4 :: mat_shape(1:2)
    
            ! Todo: trace only available for square matrices
            mat_shape = shape(entries)
            NN = mat_shape(1)
            trace = (0.d0,0.d0)
            do ii = 1, NN
                trace = trace + entries(ii,ii)
            end do
        end function mat_trace
    
        function mat_adj(entries) result(adj_mat)
            implicit none
            complex*16, intent(in) :: entries(:,:)
            complex*16, allocatable :: adj_entries(:,:)
            type(complex_matrix) :: adj_mat
            integer*4 :: ii, jj, NN, MM
            integer*4 :: mat_shape(1:2)
    
            mat_shape = shape(entries)
            NN = mat_shape(1)
            MM = mat_shape(2)
            allocate(adj_entries(1:MM,1:NN))
            do ii = 1, NN
                do jj = 1, MM
                    adj_entries(jj, ii) = conjg(entries(ii,jj))
                end do
            end do
            adj_mat = enter_mat(adj_entries)
        end function mat_adj
    
        subroutine write_mat(new_mat,file_name)
            implicit none
            character :: file_name
            type(complex_matrix) :: new_mat
            integer*4 :: ii, NN
    
            NN = new_mat%dimensions(1)
            open(12,file = file_name)
            do ii = 1, NN
                write(12,'(*(F0.8,SP,F0.8,"i",", "))') new_mat%entries(ii,:)
            end do
            close(12)
        end subroutine

        function random_mat(N)
            ! Generates an order N complex matrix with random entries x+iy
            ! where each x and y are uniformly distributed between -1 and 1
            implicit none
            integer :: N, ii, jj
            complex*16, dimension(N,N) :: random_mat
            real*8, dimension(N,N) :: A, B
            
                call random_number(A)
                ! rescaling into an element in (-1,1)
                A = 2*A - 1
                call random_number(B)
                ! rescaling into an element in (-1,1)
                B = 2*B - 1
            
                ! populating the matrix
                do ii = 1,N
                    do jj = 1,N
                        ! cast into a complex entry
                        random_mat(ii,jj) = complex(A(ii,jj),B(ii,jj))
                    end do
                end do
            
        end function random_mat
        
        function rand_positive(N)
            ! Generates an order N positive matrix with random entries x+iy
            implicit none
            integer :: N
            complex*16, dimension(N,N) :: rand_positive
            type(complex_matrix) :: trans_mat

                rand_positive = random_mat(N)
                ! Theorem a positive semi-definite matrix A is s. t. A = dagger(M)*M
                trans_mat = mat_adj(rand_positive)
                rand_positive = matmul(trans_mat%entries,rand_positive)
            
        end function rand_positive
        
        
end module