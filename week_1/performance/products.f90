program matrix_products

implicit none
real*8 :: cumulative
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

!------------------------------------------!
! matrix product, TODO: replace by function!
!------------------------------------------!
do ii = 1,M
	do jj = 1, P
		cumulative = 0.d0
		do kk = 1, N
			cumulative = cumulative + A(ii,kk)*B(kk,jj)
		end do
		C(ii,jj) = cumulative
	end do
end do

do ii = 1,M
	write(*,*) C(ii,:)
end do

read(*,*) ii

end program matrix_products

