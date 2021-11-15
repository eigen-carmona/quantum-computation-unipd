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

end module