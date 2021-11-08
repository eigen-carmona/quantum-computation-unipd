program matrix_products

implicit none

interface
	subroutine mat_mul_1(mat_a,mat_b,mat_c)
		real*8, intent(out):: mat_c(:,:)
		real*8, intent(in):: mat_a(:,:), mat_b(:,:)
	end subroutine
end interface

interface
	subroutine mat_mul_2(mat_a,mat_b,mat_c)
		real*8, intent(out):: mat_c(:,:)
		real*8, intent(in):: mat_a(:,:), mat_b(:,:)
	end subroutine
end interface

real*8, allocatable :: A(:,:), B(:,:), C(:,:), D(:,:), E(:,:)
integer*4 :: iterations, max_N, ii, jj
real*4 :: start, finish
real*4, allocatable :: timing(:,:)

write(*,*) "Enter the timing iterations and up to how many rows the test matrices should have"
read(*,*) iterations, max_N

! Allocate the timing array
allocate(timing(2:max_N,0:3))

! Generating two square random matrices, from 2 x 2 to max_N x max_n
do ii = 2, max_N
	! "input size"
	timing(ii,0) = ii

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
		call mat_mul_1(A,B,C)
	end do
	call cpu_time(finish)
	timing(ii,1) = (finish-start)/iterations

	! Measuring performance of the second product definition
	D = transpose(B)
	E = transpose(A)
	call cpu_time(start)
	do jj = 1, iterations
		call mat_mul_1(D,E,C)
	end do
	call cpu_time(finish)
	timing(ii,2) = (finish-start)/iterations

	! Measuring performance of built-in product definition
	call cpu_time(start)
	do jj = 1, iterations
		C = matmul(A,B)
	end do
	call cpu_time(finish)
	timing(ii,3) = (finish-start)/iterations

	! Release the variables for the incoming allocation
	deallocate(A)
	deallocate(B)
	deallocate(D)
	deallocate(E)
	deallocate(C)
end do

! Write the performance measurements
do ii = 2, max_N
	write(*,*) timing(ii,0:3)
end do

end program matrix_products


subroutine mat_mul_1(mat_a,mat_b,mat_c)
	implicit none
	real*8, intent(in):: mat_a(:,:), mat_b(:,:)
	real*8, intent(out) :: mat_c(:,:)
	integer*4:: M,N,O,P,ii,jj,kk
	real*8:: cumulative
	integer*4:: shape_a(2), shape_b(2)
	
	shape_a = shape(mat_a)
	shape_b = shape(mat_b)
	M = shape_a(1)
	N = shape_a(2)
	O = shape_b(1)
	P = shape_b(2)
	
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
	

subroutine mat_mul_2(mat_a,mat_b, mat_c)
	implicit none
	real*8, intent(in):: mat_a(:,:), mat_b(:,:)
	real*8, intent(out) :: mat_c(:,:)
	integer*4:: M,N,O,P,ii,jj,kk
	real*8:: cumulative
	integer*4:: shape_a(2), shape_b(2)
	
	shape_a = shape(mat_a)
	shape_b = shape(mat_b)
	M = shape_a(1)
	N = shape_a(2)
	O = shape_b(1)
	P = shape_b(2)
	
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
	

