program test_complex_matrix
use definitions
implicit none

complex*16, allocatable :: entries(:,:)
integer*4:: ii, jj, NN
real*8 :: re, im
type(complex_matrix) :: input_mat, adj_mat

write(*,*) "Enter N"
read(*,*) NN

allocate(entries(1:NN,1:NN))

do ii = 1, NN
    do jj = 1, NN
        write(*,*) "type in entry ", ii, jj
        read(*,*) re, im
        entries(ii,jj) = complex(re,im)
    end do
end do

input_mat = enter_mat(entries)

write(*,*) "Input matrix"
do ii = 1, NN
    write(*,'(*(F0.8,SP,F0.8,"i",", "))') input_mat%entries(ii,:)
end do

adj_mat = mat_adj(input_mat%entries)

write(*,*) "Adjoint matrix"
do ii = 1, NN
    write(*,'(*(F0.8,SP,F0.8,"i",", "))') adj_mat%entries(ii,:)
end do


write(*,fmt = '(A)',advance='no') "Matrix trace: "
write(*,'(F0.8,SP,F0.8,"i")') input_mat%trace
write(*,*) "Matrix dimensions: ", input_mat%dimensions

call write_mat(input_mat, "matrix.dat")

end program test_complex_matrix