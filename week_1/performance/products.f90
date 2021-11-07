program matrix_products

implicit none

interface
	function mat_mul_1(mat_a,mat_b)
		real*8, allocatable:: mat_mul_1(:,:)
		real*8, allocatable, intent(in):: mat_a(:,:), mat_b(:,:)
	end function
end interface

real*8, allocatable :: A(:,:), B(:,:), C(:,:)
integer*4 :: M, N, O, P, ii, jj, kk

write(*,*) "How many rows and columns does the leftmost matrix have?"
read(*,*) M, N

write(*,*) "How many rows and columns does the rightmost matrix have?"
read(*,*) O, P

if (N .ne. O) stop "***The number of rows in the rightmost matrix must be equal to the number of columns in the leftmost matrix.***"
! TODO: validate that the specified integers are higher than zero.

allocate(A(1:M,1:N))
allocate(B(1:O,1:P))
allocate(C(1:M,1:P))

do ii = 1,M
	write(*,*) "Specify the ",ii," row of matrix A"
	read(*,*) A(ii,:)
end do

do ii = 1,O
	write(*,*) "Specify the ",ii," row of matrix B"
	read(*,*) B(ii,:)
end do

C = mat_mul_1(A,B)

do ii = 1,M
	write(*,*) C(ii,1:P)
end do

end program matrix_products

!------------------------------------------!
! matrix product, TODO: replace by function!
!------------------------------------------!
function mat_mul_1(mat_a,mat_b)
implicit none
real*8, allocatable, intent(in):: mat_a(:,:), mat_b(:,:)
real*8, allocatable :: mat_mul_1(:,:)
integer*4:: M,N,O,P,ii,jj,kk
real*8:: cumulative
integer*4:: shape_a(2), shape_b(2)

shape_a = shape(mat_a)
shape_b = shape(mat_b)
M = shape_a(1)
N = shape_a(2)
O = shape_b(1)
P = shape_b(2)

allocate(mat_mul_1(1:M,1:P))

do ii = 1,M
	do jj = 1, P
		cumulative = 0.d0
		do kk = 1, N
			cumulative = cumulative + mat_a(ii,kk)*mat_b(kk,jj)
		end do
		mat_mul_1(ii,jj) = cumulative
	end do
end do

end function
